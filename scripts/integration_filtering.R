library(GenomicRanges)
library(stringdist)
library(igraph)
library(tidyverse)
library(spgs)

# ========================
# Functions

rc <- function(x) {toupper(spgs::reverseComplement(x))}

read_barcodes <- function(path, sample) {
  bc <- read_tsv(path, col_names = c("chr", "position", "strand", "depth_total", "barcode", "depth_barcode")) %>% 
    mutate(barcode = rc(barcode)) %>% distinct() %>% mutate(strand = ifelse(strand == "+", "-", "+"), sample = sample)
}

mode_value <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

merge_barcodes <- function(df) {
  
  # Process barcodes by length so that Hamming distances make sense.
  df <- df %>% mutate(len = nchar(barcode))
  
  merged_list <- lapply(unique(df$len), function(l) {
    df_sub <- df %>% filter(len == l)
    if(nrow(df_sub) == 0) return(NULL)
    
    # Compute pairwise Hamming distances.
    dmat <- as.matrix(stringdist::stringdistmatrix(df_sub$barcode, df_sub$barcode, method = "hamming"))
    
    # Build an undirected graph: edge if distance <= 1.
    g <- igraph::graph_from_adjacency_matrix(dmat <= 1, mode = "undirected", diag = FALSE)
    
    # Find connected components (clusters).
    comps <- igraph::components(g)
    df_sub$cluster <- comps$membership
    
    # For each cluster, select the barcode with the highest count and sum all counts.
    df_sub %>%
      group_by(cluster) %>%
      summarise(
        barcode = barcode[which.max(count)],
        count = sum(count),
        .groups = "drop"
      ) %>%
      dplyr::select(barcode, count)
  })
  
  merged_df <- bind_rows(merged_list)
  return(merged_df)
}

consolidate_barcodes_one_mismatch <- function(df) {
  barcodes <- df$barcode
  depth_barcode_vec <- df$depth_barcode
  n <- length(barcodes)
  
  # Compute all-pairs Hamming distances among barcodes
  dist_mat <- stringdist::stringdistmatrix(barcodes, barcodes, method = "hamming")
  
  # Create adjacency (TRUE if distance <= 1)
  adjacency <- dist_mat <= 1
  
  # Build graph from the adjacency matrix
  g <- igraph::graph_from_adjacency_matrix(adjacency, mode = "undirected", diag = FALSE)
  
  # Find connected components
  comps <- igraph::components(g)$membership
  
  # For each connected component, pick the barcode with the highest depth_barcode
  new_barcodes <- character(n)
  
  for (c_id in unique(comps)) {
    idx <- which(comps == c_id)
    # among barcodes in this component, pick the index of the one with max coverage
    best_idx <- idx[ which.max(depth_barcode_vec[idx]) ]
    # assign that "best" barcode label to all barcodes in this component
    new_barcodes[idx] <- barcodes[best_idx]
  }
  
  return(new_barcodes)
}

cluster_barcodes <- function(df) {
  barcode_locs_clustered <- df %>%
    group_by(chr, strand, barcode) %>%
    arrange(position, .by_group = TRUE) %>%
    # Define cluster_id: a new cluster starts whenever the gap > 500 bp
    mutate(cluster_id = cumsum(if_else(row_number() == 1, 1L, (position - lag(position)) > 500L))) %>%
    # Group by cluster
    group_by(chr, strand, barcode, cluster_id) %>%
    # Compute cluster-wide sums of depth_total and depth_barcode
    mutate(
      cluster_depth_total    = sum(depth_total, na.rm = TRUE),
      cluster_depth_total_barcode = sum(depth_barcode,   na.rm = TRUE)
    ) %>%
    # Keep only the row with the maximum depth_barcode in each cluster
    slice_max(depth_barcode, with_ties = FALSE) %>%
    ungroup()
  
  barcode_locs_merged <- barcode_locs_clustered %>%
    group_by(chr, position, strand, cluster_id) %>%
    # Merge barcodes that differ by ≤1 mismatch
    mutate(
      barcode = consolidate_barcodes_one_mismatch(cur_data())
    ) %>%
    ungroup() %>%
    group_by(chr, position, strand, cluster_id, barcode) %>%
    summarize(
      depth_barcode = sum(depth_barcode),
      depth_total = max(depth_total),
      .groups = "drop") %>% 
    group_by(barcode) %>% 
    mutate(fract_barcode = depth_barcode/sum(depth_barcode), n_mappings = dplyr::n()) %>%
    ungroup()
  
  write_tsv(barcode_locs_merged, paste0("./prc_data/bc/sites_", df$sample[1], ".tsv"))
}


# ========================
# Load files

chromosome_sizes <- read_tsv("./refgenome/hg38_len.tsv", col_names = c("chr", "end")) %>%
  mutate(start = 1, chr = paste0("chr", chr), chr_index = rev(1:24)) %>% dplyr::select(chr, start, end, chr_index)

centromeres <- read_tsv("./hg38_centromeres.tsv", col_names = c("chr", "start", "end")) %>% mutate(chr_index = rev(1:24))
whitelist <- read_tsv("./prc_data/whitelisted_barcodes.tsv.gz") %>% pull(barcode)

# Read in barcodes
bc_PB_50k_R1 <- read_barcodes("./prc_data/bc/50k_d0_map.out.txt", "PB_50k_R1")
bc_PB_5k_R1_1 <- read_barcodes("./prc_data/bc/5K_PBhapl_mapping.out.txt", "PB_5k_R1")
bc_PB_5k_R1_2 <- read_barcodes("./prc_data/bc/5K_PBhapl_mapping_2.out.txt", "PB_5k_R1")
min_LV_5k_P1_1 <- read_barcodes("./prc_data/bc/MIN_L5K_d0_map_P1_1.out.txt", "LV_5k_P1_1")
min_LV_5k_P1_2 <- read_barcodes("./prc_data/bc/MIN_L5K_d0_map_P1_2.out.txt", "LV_5k_P1_2")
min_LV_5k_P2_1 <- read_barcodes("./prc_data/bc/MIN_L5K_d0_map_P2_1.out.txt", "LV_5k_P2_1")
min_LV_5k_P2_2 <- read_barcodes("./prc_data/bc/MIN_L5K_d0_map_P2_2.out.txt", "LV_5k_P2_2")
GFP_LV_5k_P1_1 <- read_barcodes("./prc_data/bc/MIN_L5K_d0_map_P1_1_GFP.out.txt", "LV_GFP_5k_P1_1")
GFP_LV_5k_P1_2 <- read_barcodes("./prc_data/bc/MIN_L5K_d0_map_P1_2_GFP.out.txt", "LV_GFP_5k_P1_2")
GFP_LV_5k_P2_1 <- read_barcodes("./prc_data/bc/MIN_L5K_d0_map_P2_1_GFP.out.txt", "LV_GFP_5k_P2_1")
GFP_LV_5k_P2_2 <- read_barcodes("./prc_data/bc/MIN_L5K_d0_map_P2_2_GFP.out.txt", "LV_GFP_5k_P2_2")
bc_1kV3 <- read_barcodes("./prc_data/bc/1kV3_map_minion.out.txt", "1kV3")


# ========================
# Execute main functions

# Cluster barcodes
sites_PB_50k_raw <- cluster_barcodes(bc_PB_50k_R1)
sites_PB_5k_R1_raw <- cluster_barcodes(bind_rows(bc_PB_5k_R1_1, bc_PB_5k_R1_2))
sites_MIN_LV_5k_P1_raw <- cluster_barcodes(bind_rows(min_LV_5k_P1_1, min_LV_5k_P1_2))
sites_MIN_LV_5k_P2_raw <- cluster_barcodes(bind_rows(min_LV_5k_P2_1, min_LV_5k_P2_2))
sites_GFP_LV_5k_P1_raw <- cluster_barcodes(bind_rows(GFP_LV_5k_P1_1, GFP_LV_5k_P1_2))
sites_GFP_LV_5k_P2_raw <- cluster_barcodes(bind_rows(GFP_LV_5k_P2_1, GFP_LV_5k_P2_2))
sites_1kV3_minion_raw <- cluster_barcodes(bc_1kV3)

# filter integration sites
sites_MIN_LV_5k_P1 <- sites_MIN_LV_5k_P1_raw %>% filter(depth_barcode > 10, fract_barcode > 0.9, depth_barcode/depth_total > 0.8)
sites_MIN_LV_5k_P2 <- sites_MIN_LV_5k_P2_raw %>% filter(depth_barcode > 10, fract_barcode > 0.9, depth_barcode/depth_total > 0.8)
sites_GFP_LV_5k_P1 <- sites_GFP_LV_5k_P1_raw %>% filter(depth_barcode > 0, fract_barcode > 0.9, depth_barcode/depth_total > 0.8)
sites_GFP_LV_5k_P2 <- sites_GFP_LV_5k_P2_raw %>% filter(depth_barcode > 0, fract_barcode > 0.9, depth_barcode/depth_total > 0.8)
sites_PB_5k_R1 <- sites_PB_5k_R1_raw %>% filter(depth_barcode > 5, fract_barcode > 0.9, depth_barcode/depth_total > 0.8)
sites_PB_50k <- sites_PB_50k_raw %>% filter(depth_barcode > 5, fract_barcode > 0.9, depth_barcode/depth_total > 0.8)
sites_1kV3_minion <- sites_1kV3_minion_raw %>% filter(depth_barcode > 0, fract_barcode > 0.4, depth_barcode/depth_total > 0.8, barcode %in% whitelist) %>% mutate(barcode = rc(barcode))

# safe final set
write_tsv(sites_PB_5k_R1, "./prc_data/bc/sites_PB_5k_R1.tsv")
write_tsv(sites_PB_50k, "./prc_data/bc/sites_PB_50k.tsv")
write_tsv(bind_rows(sites_MIN_LV_5k_P1, sites_MIN_LV_5k_P2) %>% mutate(barcode = rc(barcode)), "./prc_data/bc/sites_L10K.tsv")
write_tsv(bind_rows(sites_GFP_LV_5k_P1, sites_GFP_LV_5k_P2) %>% mutate(barcode = rc(barcode)), "./prc_data/bc/sites_L10K_GFP.tsv")
write_tsv(sites_1kV3_minion, "./prc_data/bc/sites_1kV3.tsv")