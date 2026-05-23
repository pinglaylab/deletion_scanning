library(tidyverse)
library(fuzzyjoin)

# ==== Defining functions ====
rc <- function(x) {toupper(spgs::reverseComplement(x))}

merge_similar <- function(x) {
  combined <- x %>%
    group_by(barcode, chr_2, strand_2, sample) %>%
    arrange(end, .by_group = TRUE) %>%
    # Create a cluster identifier: increment when the gap > 5 nt
    mutate(cluster = cumsum(c(TRUE, diff(end) > 5))) %>%
    group_by(barcode, chr_2, strand_2, cluster, cell_line, day, category, replicate, ploidy, library, sample) %>%
    summarise(
      end = round(mean(end)),
      n_entries = n(), # number of entries in the cluster
      n_umi = sum(n_umi), # number of umis in the cluster
      .groups = "drop"
    )
}

mode_value <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

merge_similar_mn <- function(x) {
  combined <- x %>%
    group_by(barcode, start, chr_2, strand_2, sample) %>%
    arrange(end, .by_group = TRUE) %>%
    # Create a cluster identifier: increment when the gap > 5 nt
    mutate(cluster = cumsum(c(TRUE, diff(end) > 5))) %>%
    group_by(barcode, start, chr_2, strand_2, cluster, sample) %>%
    summarise(
      end = round(mean(end)),
      n_entries = n(), # number of entries in the cluster
      n_umi = sum(n_umi), # number of umis in the cluster
      read_name = paste(read_name, collapse = ";"),
      .groups = "drop"
    )
}

read_deletions <- function(files, name) {
  deletions <- files %>%
    set_names(tools::file_path_sans_ext(basename(.))) %>% 
    map_dfr(~ read_tsv(.x, col_names = c("chr_2", "start", "end", "id", "score", "strand_2")), .id = "sample") %>% 
    tidyr::extract(id, into = c("UMI", "barcode"), regex = ".*:(.*)_(.*)$") %>% # Extract UMI from the 'id' column
    group_by(barcode, chr_2, start, end, strand_2, sample) %>%
    summarise(n_umi = n_distinct(UMI)) %>%
    ungroup() %>% mutate(length = end - start, end = ifelse(strand_2 == "+", end, start)) %>%
    filter(length > 50) %>% # require at least 50 aligned base pairs
    separate(sample, into = c("cell_line", "day", "category", "replicate", "ploidy", "library", "id"), sep = '_') %>% dplyr::select(-id) %>%
    mutate(sample = paste(category, day, replicate))
  write_tsv(deletions, paste0("./prc_data/deletions_prc/", name, "_raw.tsv.gz"))
  return(deletions)
}

map_deletions <- function(path, beacon_sites, reference_chr = "reference") {
  sample = basename(path) %>% str_remove("_raw.tsv.gz")
  
  print(sample)
  
  deletions <- read_tsv(path)
  deletions_mapped <- deletions %>% 
    filter(chr_2 != reference_chr) %>%
    merge_similar() %>%
    fuzzyjoin::stringdist_left_join(
      beacon_sites %>% dplyr::select(-cell_line),
      by = "barcode",
      method = "hamming",
      max_dist = 1,
      distance_col = "barcode_dist"
    ) %>%
    # keep best match per deletion row (guards against duplicate matches)
    group_by(barcode.x, chr_2, strand_2, cluster, cell_line, day, category, replicate, ploidy, library, sample, end, n_entries, n_umi) %>%
    slice_min(order_by = barcode_dist, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    filter(chr == chr_2) %>%
    mutate(barcode = barcode.y) %>%
    merge_similar() %>%
    left_join(beacon_sites, by = c("barcode", "cell_line")) %>%
    mutate(start = position, length = abs(end - start)) %>%
    dplyr::select(-position) %>%
    filter(length > 100, length < 500000, (strand == "+" & start > end  & strand_2 == "+") | (strand == "-" & end > start  & strand_2 == "-"))
  write_tsv(deletions_mapped, paste0("./prc_data/deletions_prc/", sample, "_mapped.tsv.gz"))
}

find_unique_deletions <- function(deletions_mapped, gaplength = 300) {
  deletions_unique <- deletions_mapped %>%
    # 1) group constant columns
    group_by(barcode, depth_barcode, chr, start, strand, chr_2, sample, category, day, ploidy, replicate, cell_line) %>%
    # 2) Make clusters whenever a gap exceedes 500
    arrange(end, .by_group = TRUE) %>%
    mutate(cluster = cumsum((end - lag(end, default = dplyr::first(end))) > gaplength)) %>%
    ungroup() %>%
    # 3_ Collapse each cluster
    group_by(barcode, depth_barcode, chr, start, strand, chr_2, sample, category, day, ploidy, replicate, cell_line, cluster) %>%
    summarise(end = round(median(end)), length = round(median(length)), n_umi = sum(n_umi), .groups  = "drop")
}

# ==== Loading files ====

sites_PB_5k_R1 <- read_tsv("./prc_data/bc/sites_PB_5k_R1.tsv") %>% mutate(cell_line = "5K")
sites_PB_min5k <- read_tsv("./prc_data/bc/sites_PB_min5k.tsv") %>% mutate(cell_line = "min5k")
sites_L10K <- read_tsv("./prc_data/bc/sites_L10K.tsv") %>% mutate(cell_line = "L10K")
sites_L10K_GFP <- read_tsv("./prc_data/bc/sites_L10K_GFP.tsv") %>% mutate(cell_line = "L10K")
sites_1KV3 <- read_tsv("./prc_data/bc/sites_1KV3.tsv") %>% mutate(cell_line = "1KV3")

bed_50k <- list.files(path = "./prc_data/deletions", pattern = "^50k.*\\.bed.gz$", full.names = TRUE)
bed_5k <- list.files(path = "./prc_data/deletions", pattern = "^5K_d15.*\\.bed.gz$", full.names = TRUE)
bed_5k_d3 <- list.files(path = "./prc_data/deletions", pattern = "^5K_d3.*\\.bed.gz$", full.names = TRUE)
bed_min5k <- list.files(path = "./prc_data/deletions", pattern = "^min5k_d.*\\.bed.gz$", full.names = TRUE)
bed_1KV3 <- list.files(path = "./prc_data/deletions", pattern = "^1KV3.*\\.bed.gz$", full.names = TRUE)
bed_1KV3_haploid <- list.files(path = "./prc_data/deletions", pattern = "^1KV3_d18.*\\.bed.gz$", full.names = TRUE)
bed_c4 <- list.files(path = "./prc_data/deletions", pattern = "^c4_GFP.*\\.bed.gz$", full.names = TRUE)

bed_L10k_GFP <- list.files(path = "./prc_data/deletions", pattern = "^L10K.*_GFP\\.bed.gz$", full.names = TRUE)
bed_L10k_d3 <- list.files(path = "./prc_data/deletions", pattern = "^L10K_d3.*_target\\.bed.gz$", full.names = TRUE)
bed_L10k_d15 <- list.files(path = "./prc_data/deletions", pattern = "^L10K_d15.*_target\\.bed.gz$", full.names = TRUE)
bed_L10k_d25 <- list.files(path = "./prc_data/deletions", pattern = "^L10K_d25.*_target\\.bed.gz$", full.names = TRUE)
bed_L10k_d25_GFP <- list.files(path = "./prc_data/deletions", pattern = "^L10K_d25.*_GFP\\.bed.gz$", full.names = TRUE)

# ==== Identify deletions ====
deletions_1KV3_d4 <- read_deletions(bed_1KV3, "deletions_1KV3")
deletions_1KV3_d18 <- read_deletions(bed_1KV3_haploid, "deletions_1KV3_d18_haploid")
deletions_L10K_d3 <- read_deletions(bed_L10k_d3, "deletions_L10k_d3")
deletions_L10K_d15 <- read_deletions(bed_L10k_d15, "deletions_L10k_d15")
deletions_L10K_GFP <- read_deletions(bed_L10k_GFP, "deletions_L10k_GFP")
deletions_50k <- read_deletions(bed_50k, "deletions_50k")
deletions_5k_d3 <- read_deletions(bed_5k_d3, "deletions_5k_d3") %>% mutate(sample = paste(ploidy, day, replicate))
deletions_5k_d15 <- read_deletions(bed_5k, "deletions_5k_d15") %>% mutate(sample = paste(ploidy, day, replicate))

# ==== Map deletions to beacon barcodes ====

# Map deletions to beacon barcodes
deletions_50k_mapped <- map_deletions("./prc_data/deletions_prc/deletions_50k_raw.tsv.gz", sites_PB_50k)
deletions_5k_mapped <- map_deletions("./prc_data/deletions_prc/deletions_5k_raw.tsv.gz", sites_PB_5k_R1)
deletions_L10K_mapped <- map_deletions("./prc_data/deletions_prc/deletions_L10K_raw.tsv.gz", sites_L10K, reference_chr = "lenti_genome")
deletions_L10K_GFP_mapped <- map_deletions("./prc_data/deletions_prc/deletions_L10K_GFP_raw.tsv.gz", sites_L10K_GFP, reference_chr = "LV_egfp_full")
deletions_1KV3_mapped <- map_deletions("./prc_data/deletions_prc/deletions_1KV3_raw.tsv.gz", sites_1KV3)
deletions_1KV3_late_mapped <- map_deletions("./prc_data/deletions_prc/deletions_1KV3_d18_haploid_raw.tsv.gz", sites_1KV3)
deletions_min5k_mapped <- map_deletions("./prc_data/deletions_prc/deletions_min5k_raw.tsv.gz", sites_PB_min5k)
deletions_c4_mapped <- map_deletions("./prc_data/deletions_prc/deletions_c4_raw.tsv.gz", site_c4)

# ==== Filtering for unique deletions ====

deletions_1KV3_unique <- find_unique_deletions(deletions_1KV3_mapped)
deletions_1KV3_late_unique <- find_unique_deletions(deletions_1KV3_late_mapped)
deletions_L10K_unique <- find_unique_deletions(deletions_L10K_mapped)
deletions_L10K_GFP_unique <- find_unique_deletions(deletions_L10K_GFP_mapped)
deletions_5k_unique <- find_unique_deletions(deletions_5k_mapped)
deletions_min5k_unique <- find_unique_deletions(deletions_min5k_mapped)
deletions_50k_unique <- find_unique_deletions(bind_rows(deletions_50k_mapped, deletions_mn_mapped))
deletions_c4_unique <- find_unique_deletions(deletions_c4_mapped)


# ==== Deletions from minion sequencing ====

min_files <- list.files(path = "./prc_data/deletions/minion", pattern = "_deletions.tsv$", full.names = TRUE)

deletions_mn <- min_files %>%
  set_names(tools::file_path_sans_ext(basename(.))) %>% 
  map_dfr(~ read_tsv(.x), .id = "sample") %>% 
  separate(sample, into = c("ID", "day", "replicate", "category", "ploidy")) %>% 
  mutate(sample = paste(category, day, replicate), cell_line = "50k", ploidy = "diploid", library = "1") %>%
  filter(start > 178, start < 1200) %>%
  group_by(barcode, chr_2, start, end, strand_2, day, category, replicate, library, sample, cell_line, ploidy) %>%
  summarise(n_umi = n(), read_name = paste(read_name, collapse = ";")) %>%
  ungroup()
write_tsv(mutate(deletions_mn, strand_2 = ifelse(strand_2 == "+", "-", "+")), "./prc_data/deletions_prc/deletions_mn_raw.tsv.gz")

deletions_mn_mapped <- map_deletions("./prc_data/deletions_prc/deletions_mn_raw.tsv.gz", sites_PB_50k)

map_mn <- deletions_mn %>%
  group_by(barcode, chr_2, end, strand_2, sample) %>%
  filter(start == 1281) %>%
  summarise(n_umi = n(), .groups = "drop") %>%
  merge_similar() %>%
  group_by(barcode, chr_2, strand_2) %>%
  summarise(end = mode_value(end), depth_barcode = sum(n_umi)*10, depth_total = sum(n_umi)*10) %>%
  mutate(strand = ifelse(strand_2 == "+", "-", "+"), sample = "50k_d0") %>%
  dplyr::select("barcode", "chr" = "chr_2", "position" = "end", "strand", "depth_barcode", "depth_total")
write_tsv(map_mn, "./prc_data/bc/bc_minion.tsv")

deletions_mn_unique <- deletions_mn_mapped %>%
  # 1) group constant columns
  group_by(barcode, depth_barcode, chr, start, strand, chr_2, sample, category, day, ploidy, replicate, cell_line) %>%
  # 2) Make clusters whenever a gap exceedes 200
  arrange(end, .by_group = TRUE) %>%
  mutate(cluster = cumsum((end - lag(end, default = dplyr::first(end))) > 200)) %>%
  ungroup() %>%
  # 3_ Collapse each cluster
  group_by(barcode, depth_barcode, chr, start, strand, chr_2, sample, category, day, ploidy, replicate, cell_line, cluster) %>%
  summarise(end = round(median(end)), length = round(median(length)), n_umi = sum(n_umi), .groups  = "drop")

# ==== Combining and annotating deletions ====
annotation <- tibble(
  sample = c("1KV3 pre selection R1", "1KV3 pre selection R2", "L10K pre selection R1", "L10K pre selection R2", "5K pre selection R1", "5K pre selection R2",
             "1KV3 post selection haploid R1", "1KV3 post selection haploid R2", "L10K post selection haploid R1", "L10K post selection haploid R2", "5K post selection haploid R1", "5K post selection haploid R2",
             "50k pre selection R1", "50k pre selection R2", "5K pre selection R1", "5K pre selection R2",
             "50k post selection diploid R1", "50k post selection diploid R2", "5K post selection diploid R1", "5K post selection diploid R2"),
  ploidy = c(rep("haploid", 12), rep("diploid", 8)),
  selection = c(rep("pre selection", 6), rep("post selection", 6), rep("pre selection", 4), rep("post selection", 4)))

preselection <- c("L10K Cas3 d3 R1", "L10K Cas3 d3 R2", "5K mixed d3 R1", "5K mixed d3 R2", "1KV3 Cas3 d4 R1", "1KV3 Cas3 d4 R2", "50k Cas3 d3 R1", "50k Cas3 d3 R2")
postselection_haploid <- c("1KV3 haploid d18 R1", "1KV3 haploid d18 R2", "5K haploid d15 R1", "5K haploid d15 R2", "L10K GFPn d15 R1", "L10K GFPn d15 R2")
postselection_diploid <- c("5K diploid d15 R1", "5K diploid d15 R2", "50k Cas3 d15 R1", "50k Cas3 d15 R2", "50k Cas3 d27 R1", "50k Cas3 d27 R2")

# Further filtering the 5K experiment based on ploidy
barcode_ploidy <- deletions_5k_mapped %>% filter(n_umi > 5) %>% group_by(barcode, ploidy) %>% summarise(n = n()) %>% pivot_wider(names_from = ploidy, values_from = n, values_fill = 0) %>%
  mutate(sum = haploid+diploid, fraction_haploid = haploid/sum)
haploid_barcodes <- barcode_ploidy %>% filter(fraction_haploid > 0.75, sum > 4) %>% pull(barcode) %>% unique()
diploid_barcodes <- barcode_ploidy %>% filter(fraction_haploid < 0.25, sum > 4) %>% pull(barcode) %>% unique()

deletions_5k_unique_filtered <- deletions_5k_unique %>% 
  mutate(ploidy = case_when(barcode %in% haploid_barcodes ~ "haploid", barcode %in% diploid_barcodes ~ "diploid", .default = "unknown"),
         sample = paste(ploidy, day, replicate)) %>%
  filter(ploidy != "unknown")

# Now combing everything
deletions_unique <- bind_rows(deletions_L10K_unique, deletions_L10K_GFP_unique, mutate(deletions_5k_unique_filtered, sample = ifelse(day == "d3", paste("mixed", day, replicate), sample)), deletions_50k_unique, 
                              deletions_1KV3_unique, mutate(deletions_1KV3_late_unique, sample = paste(ploidy, day, replicate))) %>%
  mutate(sample = paste(cell_line, sample))

deletions_unique_selection <- deletions_unique %>% mutate(selection = case_when(
  sample %in% preselection ~ "pre selection",
  sample %in% postselection_diploid ~ "post selection diploid",
  sample %in% postselection_haploid ~ "post selection haploid"
)) %>%
  filter(!is.na(selection), length > 1000) %>%
  separate(sample, into = c("cell_line", "category", "day", "replicate"), remove = FALSE) %>%
  mutate(sample = paste(cell_line, selection, replicate)) %>%
  dplyr::select(-selection, -ploidy) %>%
  left_join(annotation, by = "sample") %>%
  mutate(identity = paste(cell_line, replicate), group = paste(cell_line, selection, ploidy)) %>%
  group_by(group, barcode) %>%
  mutate(n_deletions = n()) %>%
  ungroup()
write_tsv(deletions_unique_selection, "./prc_data/deletions_prc/deletions_unique_selection.tsv")
