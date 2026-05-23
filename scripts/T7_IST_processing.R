library(tidyverse)

setwd("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/postdoc/submission/git")

# ==== Functions ====
compute_deletion_label <- function(df, type = "translocation") {
  df %>%
    mutate(
      deletion = paste0(chr, ":", start, "-", end),
      length = if (type == "translocation") abs(end - start) else abs(end - transposon_start),
      start_rounded = paste0(round(start, -6) / 1e6, "Mb"),
      length_rounded = if (type == "translocation") {
        paste0(round(length, -3) / 1e3, "kb")
      } else {
        paste0(round(length, -2) / 1e3, "kb")
      },
      deletion_label = paste0(chr, ":", start_rounded, " | ", length_rounded)
    ) %>%
    dplyr::select(-start_rounded, -length_rounded)
}

harmonize_intervals <- function(df, gaplength = 300) {
  df %>%
    group_by(chr, start, strand, barcode, chr_2) %>%
    arrange(end, .by_group = TRUE) %>%
    mutate(cluster = cumsum((end - lag(end, default = dplyr::first(end))) > gaplength)) %>%
    ungroup() %>%
    group_by(chr, start, strand, barcode, chr_2, cluster) %>%
    mutate(
      end = round(median(end)),
      length = round(median(length))
    ) %>%
    ungroup()
}

map_barcodes_to_clones <- function(df, max_dist = 1, name) {
  mapped_df <- df %>%
    fuzzyjoin::stringdist_left_join(
      clones_scShred_v2,
      by = "barcode",
      method = "hamming",
      max_dist = max_dist,
      distance_col = "dist"
    ) %>%
    filter(!is.na(clone)) %>%
    group_by(chr, position, start, chr_2, end, strand, strand_2, cell, barcode.x, clone) %>%
    slice_min(dist, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(cell, barcode = barcode.y, clone, chr, position, start, chr_2, end, strand, strand_2, n_umi)
  write_tsv(mapped_df, paste0("./scShred/deletions/", name, ".tsv.gz"))
  return(mapped_df)
}

# ==== Loading in files
passing_cells <- read_tsv("./scShred/passing_cells_v2.tsv")
clones_scShred <- read_tsv("./clones_scShred.tsv")
clones_scShred_v2 <- clones_scShred %>% mutate(barcode = paste0(substr(barcode, 1, 4), substr(barcode, 6, nchar(barcode))))

sc_deletions_v2_raw <- read_tsv("./scShred/deletions/sc_deletions.tsv.gz") %>%
  mutate(barcode = paste0(substr(barcode, 1, 4), substr(barcode, 6, nchar(barcode)))) # remove the N in positions of low diversity on the sequencing run

# ==== Create objects  ====
sc_mapped_v2 <- sc_deletions_v2 %>%
  filter(cell %in% passing_cells) %>%
  map_barcodes_to_clones(name = "scShred_v2")

barnyard_v2 <- sc_mapped_v2 %>% group_by(cell, clone) %>%
  summarise(n_umi = sum(n_umi), .groups = "drop") %>%
  group_by(cell) %>%
  mutate(purity = max(n_umi)/sum(n_umi), dominant_clone = clone[which.max(n_umi)], total_umi = sum(n_umi)) %>%
  filter(total_umi > 99, total_umi < 2000) %>%
  pivot_wider(names_from = clone, values_from = n_umi, values_fill = 0) %>%
  mutate(clone_assignment = case_when(
    dominant_clone == "Clone1" & purity > 0.95 ~ "Clone1",
    dominant_clone == "Clone2" & purity > 0.7 ~ "Clone2",
    .default = "Unassigned"))

barnyard_long_v2 <- barnyard_v2 %>%
  filter(clone_assignment != "Low UMI") %>%
  dplyr::select(cell, clone_assignment, Clone1, Clone2) %>%
  pivot_longer(
    cols = c("Clone1", "Clone2"),
    names_to = "clone",
    values_to = "value",
  )

pairwise_df_v2 <- barnyard_long_v2 %>%
  inner_join(
    barnyard_long_v2,
    by = c("cell", "clone_assignment"),
    suffix = c("_x", "_y")
  ) %>%
  filter(clone_x < clone_y)   # keeps each pair only once

passing_clones_v2 <- barnyard_v2 %>%
  filter(!clone_assignment %in% c("Unassigned", "Low UMI")) %>%
  dplyr::select(cell, "Clone" = "clone_assignment") %>%
  distinct()
write_tsv(passing_clones_v2, "./scShred/passing_clones_v2.tsv.gz")

sc_deletions_v2 <- sc_mapped_v2 %>% 
  group_by(cell) %>%
  mutate(n_umi_cell = sum(n_umi), length = abs(end - start)) %>%
  ungroup() %>%
  filter(chr_2 == chr) %>%
  mutate(start = position, length = ifelse(chr == chr_2, abs(end - start), -1)) %>%
  filter(length > 100 | length == -1, length < 500000, (strand == "-" & start > end  & strand_2 == "-") | (strand == "+" & end > start  & strand_2 == "+")) %>%
  harmonize_intervals() %>%
  group_by(barcode, chr, start, strand, chr_2, end, clone, cell) %>%
  summarise(n_umi  = sum(n_umi), n_umi_cell = first(n_umi_cell), .groups  = "drop") %>%
  left_join(passing_clones_v2, by = "cell") %>%
  mutate(assignment = ifelse(clone == Clone, "correct", "false")) %>%
  compute_deletion_label() %>%
  group_by(deletion) %>%
  mutate(n_cells = n_distinct(cell)) %>%
  ungroup()
write_tsv(sc_deletions_v2, "./scShred/scShred_v2_deletions_prc.tsv.gz")

sc_deletions_filtered_v2 <- sc_deletions_v2 %>% 
  filter(n_umi > 10, n_umi_cell < 2000, assignment == "correct", chr == chr_2) %>%
  group_by(deletion_label) %>%
  mutate(n_cells = n_distinct(cell)) %>%
  ungroup()
write_tsv(sc_deletions_filtered_v2, "./scShred/scShred_v2_deletions_filtered.tsv.gz")
