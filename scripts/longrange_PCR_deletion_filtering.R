library(igraph)

# ==== Functions ====
mode_value <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

collapse_umis_1mm <- function(umis) {
  # Return a vector of "collapsed UMI ids" same length as umis
  umis <- as.character(umis)
  reps <- character(0)
  group_id <- integer(length(umis))
  
  for (i in seq_along(umis)) {
    u <- umis[i]
    if (length(reps) == 0) {
      reps <- u
      group_id[i] <- 1L
    } else {
      d <- stringdist(u, reps, method = "hamming")
      hit <- which(d <= 1L)[1]
      if (length(hit) == 0 || is.na(hit)) {
        reps <- c(reps, u)
        group_id[i] <- length(reps)
      } else {
        group_id[i] <- hit
      }
    }
  }
  group_id
}

collapse_start_end <- function(df, start_tol = 10, end_tol = 100) {
  # Handle empty groups gracefully
  if (nrow(df) == 0L) {
    return(data.frame(start = numeric(), end = numeric(), n_umi = integer()))
  }
  
  # 1. Aggregate identical consensus coordinates to establish our seeds.
  # This counts how many UMIs support each exact start/end pair.
  df_agg <- df %>%
    group_by(start, end) %>%
    summarise(n_umi = n(), .groups = "drop") %>%
    arrange(desc(n_umi))
  
  # Handle single-coordinate groups after aggregation
  if (nrow(df_agg) == 1L) return(df_agg)
  
  n <- nrow(df_agg)
  cluster_id <- rep(NA_integer_, n)
  current_cluster <- 1L
  
  # 2. Greedy clustering based on the aggregated seeds
  for (i in seq_len(n)) {
    if (!is.na(cluster_id[i])) next # Skip if already assigned
    
    # i is our new seed
    cluster_id[i] <- current_cluster
    
    # Find unclustered candidates
    unclustered <- which(is.na(cluster_id))
    
    if (length(unclustered) > 0L) {
      # Check tolerances against the SEED's coordinates
      is_close <- abs(df_agg$start[unclustered] - df_agg$start[i]) <= start_tol &
        abs(df_agg$end[unclustered] - df_agg$end[i]) <= end_tol
      
      # Assign matches to the current cluster
      cluster_id[unclustered[is_close]] <- current_cluster
    }
    current_cluster <- current_cluster + 1L
  }
  
  df_agg$cluster <- cluster_id
  
  # 3. Final aggregation of the clusters
  df_agg %>%
    group_by(cluster) %>%
    summarise(
      # The seed is always the first item in the cluster due to the arrange() step
      start = dplyr::first(start),
      end   = dplyr::first(end),
      n_umi = sum(n_umi),
      .groups = "drop"
    ) %>%
    dplyr::select(-cluster)
}

GFPn2k_deletions <- read_tsv("./prc_data/delamp_c4_GFPn2k.tsv.gz")

GFPn2k_deletions_UMI <- GFPn2k_deletions %>%
  filter(chr_2 == "chr4", start < 1200, strand_2 == "+") %>%
  mutate(umi_group = collapse_umis_1mm(UMI))

UMI_summarised <- GFPn2k_deletions_UMI %>% group_by(umi_group) %>%
  summarise(start = mode_value(start), end = mode_value(end), chr_2 = mode_value(chr_2), .groups = "drop") 

collapsed <- UMI_summarised %>%
  group_by(chr_2) %>%
  group_modify(~ collapse_start_end(.x, start_tol = 20, end_tol = 100)) %>%
  ungroup()

GFPn2k_deletions_filtered2 <- collapsed %>% filter(n_umi > 100) %>%
  group_by(end) %>%
  slice_max(n_umi) %>%
  group_by(start) %>%
  slice_max(n_umi) %>%
  mutate(length = 113302264-end, start = 113302264) %>%
  filter(length > 100, length < 5000000)

write_tsv(GFPn2k_deletions_filtered2, "./prc_data/deletions/delamp_c4_GFPn2k_filtered.tsv.gz")
