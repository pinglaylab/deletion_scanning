library(GenomicRanges)
library(tidyverse)
library(rtracklayer)

setwd("/Users/jonas.koeppel/Library/CloudStorage/OneDrive-Personal/postdoc/submission/git")

# Load in deletions and annotation files
deletions_final <- read_tsv("./prc_data/deletions_prc/deletions_final.tsv")
emission_annotations_HAP1 <- read_tsv("./features/ChromHMM_model/emission_annotation.txt")

# Fractional features
TAD_boundaries <- read_tsv("./features/4DNFIL8KHJD9.bed.gz", col_names = c("chr", "start", "end", "boundary_type", "boundary_strength")) %>%  mutate(segment = "TAD boundary") %>% dplyr::select(chr, start, end, segment) %>% GRanges()
constraints <- read_tsv("./features/constraint_z_genome_1kb.qc.download.txt.gz")
constraints_z3 <- filter(constraints, z > 3) %>% mutate(segment = "Mutation constraint (z > 3)") %>% dplyr::select(chrom, start, end, segment) %>% GRanges()
LADs <- read_tsv("./features/4DNFI4VAK5M4.bed.gz", col_names = c("chr", "start", "end")) %>% mutate(segment = "LAD") %>% GRanges()
chromHMM_segments_HAP1 <- read_tsv("./features/ChromHMM_model/HAP1_15_segments.bed.gz", col_names = c("chr", "start", "end", "emission")) %>% 
  left_join(emission_annotations_HAP1, by = "emission") %>% dplyr::select(chr, start, end, "segment" = "name") %>% GRanges()
AB_compartments <- rtracklayer::import("./features/4DNFIAXG3ZPL.bw") %>% as.data.frame() %>%
  mutate(segment = ifelse(score > 0, "A compartment", "B compartment")) %>% dplyr::select(seqnames, start, end, segment) %>% filter(!is.na(segment)) %>% GRanges()
replication_timing <- read_tsv("./features/GSE190117_blHAP1_wt_RT_50k_hg38_preqnorm.bedgraph.gz", col_names = c("chr", "start", "end", "reptime")) %>%
  mutate(segment = ifelse(reptime > 0, "Earlier replication", "Later replication")) %>% dplyr::select(chr, start, end, segment) %>% filter(!is.na(segment)) %>% GRanges()

# Gene level features
gtf_df  <- read_tsv("./features/gtf_df.tsv.gz") # already filtered for Appris transcripts
gene_structure_features <- filter(gtf_df, type %in% c("exon", "start_codon", "stop_codon", "five_prime_utr", "three_prime_utr"), !is.na(type)) %>% 
  transmute(seqnames, start, end, segment = str_to_sentence(str_replace_all(type, "_", " "))) %>% GRanges()
gene_biotype_features <- filter(gtf_df, type == "gene", gene_biotype != "rRNA", !is.na(gene_biotype)) %>% 
  transmute(seqnames, start, end, segment = str_replace_all(gene_biotype, "_", " ")) %>% mutate(segment = if_else(segment == "protein coding", "Protein coding", segment)) %>% GRanges()

# Pre-computed features
conservation <- read_tsv("./features/deletions_phyloP.tsv.gz", col_names = c("deletion_id", "length", "covered", "sum", "mean0", "Conservation 470 way")) %>% dplyr::select(deletion_id, `Conservation 470 way`)
gc <- read_tsv("./features/deletions_gc.tsv.gz", col_names = c("deletion_id", "length", "covered", "sum", "mean0", "GC content")) %>% dplyr::select(deletion_id, `GC content`)

annotate_with_bed <- function(deletions, states_gr) {
  # ── 1.  GRanges for deletions ───────────────────────────────────────
  del_gr <- GRanges(
    seqnames      = deletions$chr,
    ranges        = IRanges(deletions$start, deletions$end),
    deletion_id   = deletions$deletion_id)
  
  # ── 2.  pair-wise overlaps (ignoring strand) ────────────────────────
  hits <- findOverlaps(del_gr, states_gr, ignore.strand = TRUE)
  
  # intersection length for every hit
  ov_len <- width(pintersect(del_gr[queryHits(hits)],
                             states_gr[subjectHits(hits)]))
  
  # build a data frame of overlaps
  ov_tbl <- tibble(
    deletion_id = mcols(del_gr)$deletion_id[queryHits(hits)],
    state       = mcols(states_gr)$segment[subjectHits(hits)],
    bp          = ov_len
  )
  
  # ── 3.  add deletion length & turn into percentages ─────────────────
  del_len_tbl <- tibble(
    deletion_id = mcols(del_gr)$deletion_id,
    del_len     = width(del_gr)
  )
  
  annot_long <- ov_tbl %>%
    group_by(deletion_id, state) %>%
    summarise(bp = sum(bp), .groups = "drop") %>%
    left_join(del_len_tbl, by = "deletion_id") %>%
    mutate(frac = bp / del_len)
  
  # ── 4.  (optional) wide format: one column per state ────────────────
  annot_wide <- annot_long %>%
    dplyr::select(deletion_id, state, frac) %>%
    pivot_wider(names_from = state, values_from = frac, values_fill = 0)
  return(annot_wide)
}

annotate_with_bed_binary <- function(deletions, states_gr) {
  # 1. GRanges for deletions
  del_gr <- GRanges(
    seqnames    = deletions$chr,
    ranges      = IRanges(deletions$start, deletions$end),
    deletion_id = deletions$deletion_id
  )
  
  # 2. Pair-wise overlaps
  hits <- findOverlaps(del_gr, states_gr, ignore.strand = TRUE)
  
  # 3. Build overlap table and collapse to binary presence/absence
  ov_tbl <- tibble(
    deletion_id = mcols(del_gr)$deletion_id[queryHits(hits)],
    state       = mcols(states_gr)$segment[subjectHits(hits)],
    overlap     = 1L
  ) %>%
    distinct(deletion_id, state, .keep_all = TRUE)
  
  # 4. Wide format: one column per state, 1 if any overlap, else 0
  annot_wide <- ov_tbl %>%
    pivot_wider(
      names_from = state,
      values_from = overlap,
      values_fill = 0
    )
  
  # 5. Make sure deletions with no overlaps are retained
  all_deletions <- tibble(deletion_id = mcols(del_gr)$deletion_id)
  
  annot_wide <- all_deletions %>%
    left_join(annot_wide, by = "deletion_id") %>%
    mutate(across(-deletion_id, ~ tidyr::replace_na(., 0L)))
  
  return(annot_wide)
}

annot_chromHMM <- annotate_with_bed(deletions_final, chromHMM_segments_HAP1)
annot_LADs <- annotate_with_bed(deletions_final, LADs)
annot_AB <- annotate_with_bed(deletions_final, AB_compartments)
annot_reptime <- annotate_with_bed(deletions_final, replication_timing)
annot_gene_structure <- annotate_with_bed_binary(deletions_final, gene_structure_features)
annot_gene_biotype <- annotate_with_bed_binary(deletions_final, gene_biotype_features)
annot_z3 <- annotate_with_bed_binary(deletions_final, constraints_z3)
annot_TAD_boundaries <- annotate_with_bed_binary(deletions_final, TAD_boundaries)

annot_wide <- list(
  annot_chromHMM,
  annot_TAD_boundaries,
  annot_z3,
  annot_LADs,
  annot_AB,
  annot_gene_structure,
  annot_gene_biotype,
  annot_reptime) %>%
  purrr::reduce(full_join, by = "deletion_id") %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

# ── Merge back onto the original deletion table ─────
deletions_annot <- deletions_final %>%
  left_join(annot_wide, by = "deletion_id") %>%
  left_join(conservation, by = "deletion_id") %>%
  left_join(gc, by = "deletion_id") %>%
  left_join(dplyr::select(selection, barcode, "chr" = "seqnames", start, end, sample, "Coding gene dispensability (CRISPR)" = "min_lfc", "Coding gene dispensability (Depmap)" = "min_chronos"), by = c("barcode", "chr", "start", "end", "sample")) %>%
  left_join(dplyr::select(selection_lncRNA, barcode, "chr" = "seqnames", start, end, sample, "lncRNA dispensability (CRISPR)" = "min_lfc"), by = c("barcode", "chr", "start", "end", "sample")) %>%
  mutate(`SV length (log10)` = log10(length), `Coding gene dispensability (CRISPR)` = ifelse(is.na(`Coding gene dispensability (CRISPR)`), 0, `Coding gene dispensability (CRISPR)`), `lncRNA dispensability (CRISPR)` = ifelse(is.na(`lncRNA dispensability (CRISPR)`), 0, `lncRNA dispensability (CRISPR)`)) %>%
  dplyr::select(-length)

write_tsv(deletions_annot, "./prc_data/deletions/annotated_deletions.tsv.gz")

