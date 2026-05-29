library(tidyverse)
library(Rsamtools)
library(GenomicAlignments)
library(Biostrings)
library(fuzzyjoin)

setwd("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/postdoc/submission/git")

chr_list <- sprintf("chr%s",c(seq(1,22,1), "X", "Y"))

# ==== Define functions ====
rc <- function(x) {toupper(spgs::reverseComplement(x))}

mode_value <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Function to extract soft-clipped sequence from a CIGAR string and read sequence.
extractSoftClip <- function(cigar, seq) {
  seq <- as.character(seq)
  # Extract CIGAR operations (letters) and their lengths (numbers)
  ops <- unlist(regmatches(cigar, gregexpr("[A-Z]", cigar)))
  lens <- as.numeric(unlist(regmatches(cigar, gregexpr("[0-9]+", cigar))))
  
  left_clip <- ""
  right_clip <- ""
  
  # If the first operation is soft clipping ("S"), extract from the start.
  if (length(ops) > 0 && ops[1] == "S") {
    left_clip <- substr(seq, 1, lens[1])
  }
  # If the last operation is soft clipping ("S"), extract from the end.
  if (length(ops) > 0 && tail(ops, 1) == "S") {
    total_len <- nchar(seq)
    right_clip <- substr(seq, total_len - tail(lens, 1) + 1, total_len)
  }
  
  # Combine left and right soft clips with a dash if both exist.
  if (left_clip != "" && right_clip != "") {
    soft_clip <- paste0(left_clip, "-", right_clip)
  } else if (left_clip != "") {
    soft_clip <- left_clip
  } else if (right_clip != "") {
    soft_clip <- right_clip
  } else {
    soft_clip <- ""
  }
  
  return(soft_clip)
}

# Function to create a data frame from the scanBam output
process_bam <- function(path, sample) {
  
  # Open the BAM file and extract relevant fields.
  bam_data <- scanBam(path, param = ScanBamParam(what = c("qname", "flag", "seq", "cigar", "rname", "pos")))[[1]]
  
  # Compute alignment width on the _reference_ and thus the end coordinate
  ref_width <- GenomicAlignments::cigarWidthAlongReferenceSpace(bam_data$cigar)
  aln_end   <- bam_data$pos + ref_width - 1
  
  # Build bam tibble
  df <- tibble(
    read_name = bam_data$qname,
    flag       = bam_data$flag,
    seq        = as.character(bam_data$seq),
    cigar      = bam_data$cigar,
    
    # New columns:
    aln_chr    = as.character(bam_data$rname),
    aln_start  = bam_data$pos,
    aln_end    = aln_end
  ) %>%
    mutate(is_supplemental = bitwAnd(flag, 0x800) != 0)
  
  # Separate primary and supplemental alignments
  primary_df <- df %>% filter(!is_supplemental) %>% mutate(soft_clipped = mapply(extractSoftClip, cigar, seq))
  supplemental_df <- df %>% filter(is_supplemental) %>% mutate(supplemental = seq)
  
  # If there are multiple supplemental alignments per read_name, take the first one.
  supplemental_unique <- supplemental_df %>%
    # bring in the primary alignment’s chr & rev‐flag
    group_by(read_name) %>%
    # keep only reads where all supplemental alignments are on same chr
    filter(n_distinct(aln_chr) == 1, n_distinct(flag) == 1) %>%
    mutate(n = n()) %>%
    dplyr::select(read_name, supplemental, "sup_chr" = "aln_chr", "sup_start" = "aln_start", "sup_end" = "aln_end", n) %>%
    ungroup()
  
  # Merge primary and supplemental info by read_name. 
  final_df <- primary_df %>% 
    dplyr::select(read_name, soft_clipped, aln_chr, aln_start, aln_end) %>%
    left_join(supplemental_unique, by = "read_name") %>%
    mutate(sample = sample) %>%
    distinct() %>%
    filter(aln_chr %in% c(chr_list, "reference", "lenti_genome", "LV_egfp_full"))
  
  write_tsv(final_df, paste0("./prc_data/deletions_prc/bam_reads/bam_reads_", sample, ".tsv.gz"))
  return(final_df)
  
}

# Map deletions candidates to beacon barcodes
map_deletions_split <- function(deletions, filename, beacon_sites) {
  deletions_mapped <- deletions %>% 
    fuzzyjoin::stringdist_left_join(
      beacon_sites %>% dplyr::select(-cell_line),
      by = "barcode",
      method = "hamming",
      max_dist = 1
    ) %>%
    # keep best match per deletion row (guards against duplicate matches)
    filter(chr == chr_2) %>%
    mutate(barcode = barcode.y, transposon_start = start, start = position, length = abs(end - position)) %>%
    dplyr::select(-position, -barcode.x, -barcode.y) %>%
    filter(length > 100, length < 500000, (strand == "+" & start > end  & strand_2 == "+") | (strand == "-" & end > start  & strand_2 == "-"))
  write_tsv(deletions_mapped, paste0("./prc_data/deletions_prc/", filename, "_split_mapped.tsv.gz"))
}

# Function to only retain reads that have an associated deletion
filter_breakpoints <- function(breakpoints, deletion_expanded, end) {
  breakpoints_filtered <- breakpoints %>% 
    filter(sup_end < end) %>%
    inner_join(deletion_expanded, by = "read_name") %>%
    # Take the correct sequence based on the orientation of the cassette
    mutate(soft_clipped = ifelse(strand == "+", rc(soft_clipped), soft_clipped), soft_clipped = str_remove(soft_clipped, "-.*$")) %>% 
    separate(sample, into = c("cell_line", "day", "category", "replicate", "ploidy", "library", "target")) %>%
    mutate(sample = paste(cell_line, category, day, target, replicate))
  return(breakpoints_filtered)
}

# Function to annotate junctions from split reads
generate_junctions <- function(breakpoints) {
  junctions <- breakpoints %>% 
    mutate(length_clipped = nchar(soft_clipped), length_supplemental = nchar(supplemental), gap = length_clipped-length_supplemental,
           leftover = str_replace(soft_clipped, supplemental, "_SEC_"), overlap = str_replace(supplemental, soft_clipped, "_SEC_"),
           breakend_seq = ifelse(gap < 0, overlap, leftover)) %>%
    filter(length_clipped > 20) %>%
    separate(breakend_seq, into = c("left", "sec", "right"), sep = "_", remove = F) 
}

# Function to merge nearby junctions
collapse_junctions <- function(junctions, construct) {
  junctions_collapsed <- junctions %>% 
    group_by(barcode, start, strand, chr_2, strand_2, sample, target) %>%
    arrange(end, .by_group = TRUE) %>%
    mutate(cluster = cumsum((end - lag(end, default = dplyr::first(end))) > 200)) %>% # Make clusters that combine everything within 200 bp
    ungroup() %>%
    group_by(barcode, start, strand, chr_2, strand_2, sample, target, cluster) %>%
    summarise(end = round(median(end)), length = round(median(length)), n = n(), ins_seq = mode_value(right), length_clipped = mode_value(length_clipped), length_supplemental = mode_value(length_supplemental), 
              gap_frac = mean(gap == mode_value(gap), na.rm = T), gap = mode_value(gap), 
              transposon_start = mode_value(transposon_start), .groups = "drop") %>%
    filter(gap_frac > 0.8) %>%
    mutate(end_type = case_when(
      gap == 0             ~ "Precise",
      gap <  0             ~ "Microhomology",
      gap >  0             ~ "Insertion"),
      guide_position = case_when(construct == "lenti_minimal" ~ transposon_start - 2466, 
                                 construct == "lenti_GFP" ~ transposon_start - 2722,
                                 construct == "PB_GFP" ~ transposon_start - 573,
                                 construct == "PB_GFP_v3" ~ transposon_start - 586),
      resection = ifelse(guide_position < 4, "resected", "non_resected"))
  return(junctions_collapsed)
}

# ==== Loading in files ====
sites_PB_5k_R1 <- read_tsv("./prc_data/bc/sites_PB_5k_R1.tsv") %>% mutate(cell_line = "5K")
sites_PB_50k <- read_tsv("./prc_data/bc/sites_PB_50k.tsv") %>% mutate(cell_line = "50k")
sites_L10K <- read_tsv("./prc_data/bc/sites_L10K.tsv") %>% mutate(cell_line = "L10K")
sites_L10K_GFP <- read_tsv("./prc_data/bc/sites_L10K_GFP.tsv") %>% mutate(cell_line = "L10K")
sites_1KV3 <- read_tsv("./prc_data/bc/sites_1KV3.tsv") %>% mutate(cell_line = "1KV3")

# Deletion files
files_mn_split <- list.files(path = "./nanopore/raw_data/deletions", pattern = "_deletions.tsv$", full.names = TRUE)
files_L10k_split <- list.files(path = "./prc_data/deletions/", pattern = "^L10K_d.*target_deletions.tsv$", full.names = TRUE)
files_L10k_split_GFP <- list.files(path = "./prc_data/deletions/", pattern = "^L10K_d.*GFP_deletions.tsv$", full.names = TRUE)
files_50k_split <- list.files(path = "./prc_data/deletions/", pattern = "^50k_d.*_deletions.tsv$", full.names = TRUE)
files_5k_split <- list.files(path = "./prc_data/deletions/", pattern = "^5K_d.*_deletions.tsv$", full.names = TRUE)
files_1KV3_split <- list.files(path = "./prc_data/deletions/", pattern = "^1KV3.*_deletions.tsv$", full.names = TRUE)

# Bam files
files_bam_MIN <- list.files(path = "./nanopore/raw_data/deletions/", pattern = "^MIN_d.*_del.bam$", full.names = TRUE)
files_bam_L10K <- list.files(path = "./prc_data/bam/", pattern = "^L10K_d.*target_del.bam$", full.names = TRUE)
files_bam_GFP <- list.files(path = "./prc_data/bam/", pattern = "^L10K_d.*GFP_del.bam$", full.names = TRUE)
files_bam_50k <- list.files(path = "./prc_data/bam/", pattern = "^50k_d.*_del.bam$", full.names = TRUE)
files_bam_5k <- list.files(path = "./prc_data/bam/", pattern = "^5K_d.*_del.bam$", full.names = TRUE)
files_bam_1KV3 <- list.files(path = "./prc_data/bam/", pattern = "^1KV3.*_del.bam$", full.names = TRUE)

# Reading in bam files
breakpoints_1KV3 <- map_dfr(files_bam_1KV3[str_detect(files_bam_1KV3, "d4")], ~ process_bam(.x, sample = str_remove(tools::file_path_sans_ext(basename(.)), "_del")))
breakpoints_5k <- map_dfr(files_bam_5k, ~ process_bam(.x, sample = str_remove(tools::file_path_sans_ext(basename(.)), "_del")))
breakpoints_50k <- map_dfr(files_bam_50k, ~ process_bam(.x, sample = str_remove(tools::file_path_sans_ext(basename(.)), "_del")))
breakpoints_L10K_target <- map_dfr(files_bam_L10K, ~ process_bam(.x, sample = str_remove(tools::file_path_sans_ext(basename(.)), "_del")))
breakpoints_L10K_GFP <- map_dfr(files_bam_GFP, ~ process_bam(.x, sample = str_remove(tools::file_path_sans_ext(basename(.)), "_del")))
breakpoints_mn <- map_dfr(files_bam_MIN, ~ process_bam(.x, sample = str_remove(tools::file_path_sans_ext(basename(.)), "_del")))

# Reading in deletions
deletions_mn_split <- files_mn_split %>%
  set_names(tools::file_path_sans_ext(basename(.))) %>% 
  map_dfr(~ read_tsv(.x), .id = "sample") %>% 
  separate(sample, into = c("cell_line", "day", "replicate", "category", "library", "del")) %>% 
  mutate(sample = paste(category, day, replicate), cell_line = "50k", ploidy = "diploid") %>%
  filter(start > 178, start < 1200) %>%
  group_by(barcode, chr_2, start, end, strand_2, day, category, replicate, library, sample, cell_line, ploidy) %>%
  summarise(n_umi = n(), read_name = paste(read_name, collapse = ";"), .groups = "drop")

deletions_1KV3_split <- files_1KV3_split[str_detect(files_1KV3_split, "d4")] %>%
  set_names(tools::file_path_sans_ext(basename(.))) %>% 
  map_dfr(~ read_tsv(.x), .id = "sample") %>% 
  separate(sample, into = c("cell_line", "day", "category", "replicate", "ploidy", "library", "del")) %>% 
  mutate(sample = paste(category, day, replicate)) %>%
  filter(start > 178, start < 1200) %>%
  group_by(barcode, chr_2, start, end, strand_2, day, category, replicate, sample, cell_line) %>%
  summarise(n_umi = n(), read_name = paste(read_name, collapse = ";"), .groups = "drop")

deletions_5k_split <- files_5k_split %>%
  set_names(tools::file_path_sans_ext(basename(.))) %>% 
  map_dfr(~ read_tsv(.x), .id = "sample") %>% 
  separate(sample, into = c("cell_line", "day", "category", "replicate", "ploidy", "library", "del")) %>% 
  mutate(sample = paste(category, day, replicate)) %>%
  filter(start > 178, start < 1200) %>%
  group_by(barcode, chr_2, start, end, strand_2, day, category, replicate, library, sample, cell_line, ploidy) %>%
  summarise(n_umi = n(), read_name = paste(read_name, collapse = ";"), .groups = "drop")

deletions_50k_split <- files_50k_split %>%
  set_names(tools::file_path_sans_ext(basename(.))) %>% 
  map_dfr(~ read_tsv(.x), .id = "sample") %>% 
  separate(sample, into = c("cell_line", "day", "category", "replicate", "ploidy", "library", "del")) %>% 
  mutate(sample = paste(category, day, replicate)) %>%
  filter(start > 178, start < 1200) %>%
  group_by(barcode, chr_2, start, end, strand_2, day, category, replicate, library, sample, cell_line, ploidy) %>%
  summarise(n_umi = n(), read_name = paste(read_name, collapse = ";"), .groups = "drop")

deletions_L10K_split <- files_L10k_split %>%
  set_names(tools::file_path_sans_ext(basename(.))) %>% 
  map_dfr(~ read_tsv(.x), .id = "sample") %>% 
  separate(sample, into = c("cell_line", "day", "category", "replicate", "ploidy", "library", "target")) %>% 
  mutate(sample = paste(category, day, replicate)) %>%
  filter(start > 2263, start < 2604) %>% # alignment between T7 promoter and LTR
  group_by(barcode, chr_2, start, end, strand_2, day, category, replicate, library, sample, cell_line, ploidy) %>%
  summarise(n_umi = n(), read_name = paste(read_name, collapse = ";"), .groups = "drop")

deletions_L10K_split_GFP <- files_L10k_split_GFP %>%
  set_names(tools::file_path_sans_ext(basename(.))) %>% 
  map_dfr(~ read_tsv(.x), .id = "sample") %>% 
  separate(sample, into = c("cell_line", "day", "category", "replicate", "ploidy", "library", "del")) %>% 
  mutate(sample = paste(category, day, replicate)) %>%
  filter(start > 2241, start < 4053) %>% # alignment between T7 promoter and LTR
  group_by(barcode, chr_2, start, end, strand_2, day, category, replicate, library, sample, cell_line, ploidy) %>%
  summarise(n_umi = n(), read_name = paste(read_name, collapse = ";"), .groups = "drop")

# Map deletions to barcodes
deletions_1KV3_split_mapped <- map_deletions_split(filter(deletions_1KV3_split, day == "d4"), "deletions_1KV3_pre", sites_1KV3)
deletions_5k_split_mapped <- map_deletions_split(filter(deletions_5k_split, day == "d3"), "deletions_5k_pre", sites_PB_5k_R1)
deletions_50k_split_mapped <- map_deletions_split(filter(deletions_50k_split, day == "d3"), "deletions_50k_pre", sites_PB_50k)
deletions_L10K_split_mapped <- map_deletions_split(filter(deletions_L10K_split, day == "d3"), "deletions_L10K_pre", sites_L10K)
deletions_L10K_split_GFP_mapped <- map_deletions_split(filter(deletions_L10K_split_GFP, day == "d3"), "deletions_L10K_GFP_pre", sites_L10K_GFP)
deletions_mn_split_mapped <- deletions_mn_split %>% filter(str_detect(sample, "d3")) %>%
  stringdist_left_join(dplyr::select(sites_PB_50k, -cell_line), by = c("barcode"), method = "hamming", max_dist = 1) %>% 
  filter(chr == chr_2) %>%
  mutate(barcode = barcode.y, transposon_start = start, start = position, length = abs(end - position)) %>%
  dplyr::select(-position, -barcode.x, -barcode.y) %>%
  filter(length > 100, length < 500000, (strand == "+" & start > end  & strand_2 == "-") | (strand == "-" & end > start  & strand_2 == "+"))

# Expand the deletion table so that each row corresponds to one read
deletions_1KV3_split_expanded <- deletions_1KV3_split_mapped %>% separate_rows(read_name, sep = ";") %>% dplyr::select(barcode, read_name, start, strand, chr_2, end, strand_2, transposon_start, length)
deletions_5k_split_expanded <- deletions_5k_split_mapped %>% separate_rows(read_name, sep = ";") %>% dplyr::select(barcode, read_name, start, strand, chr_2, end, strand_2, transposon_start, length)
deletions_50k_split_expanded <- deletions_50k_split_mapped %>% separate_rows(read_name, sep = ";") %>% dplyr::select(barcode, read_name, start, strand, chr_2, end, strand_2, transposon_start, length)
deletions_L10K_split_expanded <- deletions_L10K_split_mapped %>% separate_rows(read_name, sep = ";") %>% dplyr::select(barcode, read_name, start, strand, chr_2, end, strand_2, transposon_start, length)
deletions_L10K_split_GFP_expanded <- deletions_L10K_split_GFP_mapped %>% separate_rows(read_name, sep = ";") %>% dplyr::select(barcode, read_name, start, strand, chr_2, end, strand_2, transposon_start, length)
deletions_mn_split_expanded <- deletions_mn_split_mapped %>% separate_rows(read_name, sep = ";") %>% dplyr::select(barcode, read_name, start, strand, chr_2, end, strand_2, transposon_start, length)
deletion_reads_mn <- deletions_mn_split_mapped %>% separate_rows(read_name, sep = ";") %>% pull(read_name)

# Filter breakpoints for deletions
breakpoints_1KV3_filtered <- filter_breakpoints(breakpoints_1KV3, deletions_1KV3_split_expanded, 1200)
breakpoints_5k_filtered <- filter_breakpoints(breakpoints_5k, deletions_5k_split_expanded, 1200)
breakpoints_50k_filtered <- filter_breakpoints(breakpoints_50k, deletions_50k_split_expanded, 1200)
breakpoints_L10K_target_filtered <- filter_breakpoints(breakpoints_L10K_target, deletions_L10K_split_expanded, 2546)
breakpoints_L10K_GFP_filtered <- filter_breakpoints(breakpoints_L10K_GFP, deletions_L10K_split_GFP_expanded, 3000)
breakpoints_mn_filtered <- breakpoints_mn %>% 
  filter(read_name %in% deletion_reads_mn, n == 1, sup_chr == "reference") %>%
  left_join(deletions_mn_split_expanded, by = "read_name") %>%
  # Take the correct sequence based on the orientation of the cassette
  mutate(soft_clipped = ifelse(strand == "+", rc(soft_clipped), soft_clipped), soft_clipped = str_remove(soft_clipped, "-.*$")) %>% 
  filter(str_detect(soft_clipped, "^GCTTTAAGGCC")) %>% 
  separate(sample, into = c("cell_line", "day", "category", "replicate", "ploidy", "library", "target")) %>%
  mutate(sample = paste(cell_line, category, day, target, replicate))

# Generate junction sequences
junction_seq_L10K_target <- generate_junctions(breakpoints_L10K_target_filtered)
junction_seq_L10K_GFP <- generate_junctions(breakpoints_L10K_GFP_filtered)
junction_seq_50k <- generate_junctions(breakpoints_50k_filtered)
junction_seq_5k <- generate_junctions(breakpoints_5k_filtered)
junction_seq_1KV3 <- generate_junctions(breakpoints_1KV3_filtered)
junction_seq_mn <- breakpoints_mn_filtered %>% 
  mutate(length_clipped = nchar(soft_clipped), length_supplemental = nchar(supplemental), gap = length_clipped-length_supplemental,
         leftover = str_replace(soft_clipped, supplemental, "_SEC_"), overlap = str_replace(supplemental, soft_clipped, "_SEC_"),
         breakend_seq = ifelse(gap < 0, overlap, leftover)) %>%
  filter(length_clipped > 20) %>%
  separate(breakend_seq, into = c("left", "sec", "right"), sep = "_", remove = F) %>%
  filter(!is.na(sec), left == "")

# Collapse junctions in proximity
junctions_collapsed_L10K <- bind_rows(collapse_junctions(junction_seq_L10K_target, "lenti_minimal"), collapse_junctions(junction_seq_L10K_GFP, "lenti_GFP"))
junctions_collapsed_50k <- collapse_junctions(junction_seq_50k, "PB_GFP")
junctions_collapsed_5k <- collapse_junctions(junction_seq_5k, "PB_GFP")
junctions_collapsed_1KV3 <- collapse_junctions(mutate(junction_seq_1KV3, target = "GFP"), "PB_GFP_v3")
junctions_collapsed_mn <- junction_seq_mn %>% 
  mutate(sample = paste(cell_line, target, replicate)) %>%
  group_by(barcode, start, strand, chr_2, strand_2, sample, target) %>%
  arrange(end, .by_group = TRUE) %>%
  mutate(cluster = cumsum((end - lag(end, default = dplyr::first(end))) > 200)) %>% # Make clusters that combine everything within 200 bp
  ungroup() %>%
  group_by(barcode, start, strand, chr_2, strand_2, sample, target, cluster) %>%
  summarise(end = round(median(end)), length = round(median(length)), n = n(), ins_seq = mode_value(right), length_clipped = mode_value(length_clipped), length_supplemental = mode_value(length_supplemental), 
            gap_frac = mean(gap == mode_value(gap), na.rm = T), gap = mode_value(gap), 
            transposon_start = mode_value(transposon_start), .groups = "drop") %>%
  filter(gap_frac > 0.9) %>%
  mutate(resection = ifelse(transposon_start < 576, "resected", "non_resected"),
         end_type = case_when(
           gap == 0             ~ "Precise",
           gap <  0             ~ "Microhomology",
           gap >  0             ~ "Insertion"),
         guide_position = transposon_start - 573)

# Combine everything
junctions <- bind_rows(mutate(junctions_collapsed_5k, sample = "PB v1 beacon"),
                       mutate(junctions_collapsed_50k, sample = "PB pilot Illumina"),
                       mutate(junctions_collapsed_mn, sample = "PB pilot Nanopore"),
                       mutate(junctions_collapsed_L10K, sample = "Lenti v2 beacon"),
                       mutate(junctions_collapsed_1KV3, sample = "PB v3 beacon"))
write_tsv(junctions, "./prc_data/breakends.tsv")
