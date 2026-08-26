# =============================================================================
# ATAC-seq + RNA-seq Integration Analysis
# PART 1: Data prep, DESeq2, four-quadrant integration (no clusterProfiler)
# =============================================================================
# Workflow:
#   Step 1: Extract nearby genes from annotated ATAC differential peaks
#           (split into closed / opened groups)
#   Step 2: DESeq2 differential expression analysis from HTseq count files
#   Step 4: Four-quadrant integration — ATAC direction vs RNA direction
#
# Step 3 (GO enrichmTreat2 on ATAC peaks/RNA up genes) and Step 5 (GO on Q1/Q2)
# are performed in the companion script (PART 2), which uses clusterProfiler.
# This script saves all objects needed by PART 2 to an .RData file.
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(org.Hs.eg.db)
  library(ggrepel)
  library(tidyr)
})

# -----------------------------------------------------------------------------
# Paths — edit these before running
# -----------------------------------------------------------------------------
ANNO   <- "/annotated_diff_peaks_FDR_th0.01.txt"
HTSEQ_DIR   <- "/htseq_count_out_final"   # <- set to your actual directory
OUTDIR      <- "."
#OUTDIR       <- "integration_Case_Control_output"
setwd(OUTDIR)

# HTseq filename -> sample metadata mapping
sample_info <- data.frame(
  filename  = c("case",contol),
  condition = c("control","control","case","case",
                "Treat1","Treat1","Treat2","Treat2"),
  replicate = c(1,2,1,2,1,2,1,2),
  stringsAsFactors = FALSE
)

# =============================================================================
# STEP 1: Extract nearby genes from annotated ATAC peaks
# =============================================================================
cat("\n=== STEP 1: Extracting nearby genes from ATAC differential peaks ===\n")

atac <- read.table(ANNO, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "")

atac <- atac %>%
  filter(!is.na(SYMBOL), SYMBOL != "", !is.na(ENSEMBL), ENSEMBL != "") %>%
  mutate(ENSEMBL_clean = sub("\\..*", "", ENSEMBL))  # strip version suffix

# Split by direction:
#   ATAC Fold = control - case  (g1=control, g2=case in DiffBind)
#   Fold > 0  -> control higher  -> accessibility DECREASED in case (closed)
#   Fold < 0  -> case higher    -> accessibility INCREASED in case (opened)
atac_closed <- atac %>% filter(Fold > 0) %>% arrange(desc(Fold))
atac_opened <- atac %>% filter(Fold < 0) %>% arrange(Fold)

# Collapse to one represTreat2ative peak per gene (most significant FDR)
genes_closed <- atac_closed %>%
  group_by(SYMBOL, ENSEMBL_clean) %>%
  summarise(max_Fold      = max(Fold),
            mean_Fold     = mean(Fold),
            n_peaks       = n(),
            min_FDR       = min(FDR),
            min_distToTSS = min(abs(distanceToTSS)),
            .groups = "drop") %>%
  arrange(desc(max_Fold))

genes_opened <- atac_opened %>%
  group_by(SYMBOL, ENSEMBL_clean) %>%
  summarise(min_Fold      = min(Fold),
            mean_Fold     = mean(Fold),
            n_peaks       = n(),
            min_FDR       = min(FDR),
            min_distToTSS = min(abs(distanceToTSS)),
            .groups = "drop") %>%
  arrange(min_Fold)

cat("Total raw ATAC peaks (all)   :", nrow(atac), "\n")
cat("Total raw ATAC peaks (closed):", nrow(atac_closed), "\n")
cat("Total raw ATAC peaks (opened):", nrow(atac_opened), "\n")
cat("Unique genes near CLOSED peaks:", nrow(genes_closed), "\n")
cat("Unique genes near OPENED peaks:", nrow(genes_opened), "\n")

write.table(genes_closed, "ATAC_closed_genes_case2control.txt",
            quote=FALSE, row.names=FALSE, sep="\t")
write.table(genes_opened, "ATAC_opened_genes_case2control.txt",
            quote=FALSE, row.names=FALSE, sep="\t")

# =============================================================================
# STEP 2: DESeq2 differential expression from HTseq counts (case vs Control)
# =============================================================================
cat("\n=== STEP 2: DESeq2 differential expression analysis ===\n")

read_htseq <- function(filepath, samplename) {
  df <- read.table(filepath, header=FALSE, sep="\t",
                   col.names=c("gene_id","count"), stringsAsFactors=FALSE)
  df <- df[!grepl("^__", df$gene_id), ]          # remove HTseq summary lines
  df$gene_id <- sub("\\..*", "", df$gene_id)      # strip Ensembl version suffix
  colnames(df)[2] <- samplename
  return(df)
}

# Use only control and case samples for this contrast
samples_use <- sample_info %>% filter(condition %in% c("control","case"))

count_list <- lapply(seq_len(nrow(samples_use)), function(i) {
  fpath <- file.path(HTSEQ_DIR, samples_use$filename[i])
  read_htseq(fpath, samples_use$filename[i])
})

# Merge into a single count matrix
count_mat <- Reduce(function(a,b) merge(a, b, by="gene_id", all=TRUE), count_list)
count_mat[is.na(count_mat)] <- 0

# Collapse duplicate gene_id rows (e.g. PAR_Y genes whose version-stripped
# Ensembl ID collides with the X-linked copy) by summing their counts
count_mat <- count_mat %>%
  group_by(gene_id) %>%
  summarise(across(everything(), sum), .groups = "drop") %>%
  as.data.frame()

rownames(count_mat) <- count_mat$gene_id
count_mat$gene_id <- NULL

# Keep only ENSG genes (exclude HTLV-1 gene rows)
count_mat_ensg <- count_mat[grepl("^ENSG", rownames(count_mat)), ]

col_data <- data.frame(
  row.names = samples_use$filename,
  condition = factor(samples_use$condition, levels=c("control","case"))
)

dds <- DESeqDataSetFromMatrix(countData = count_mat_ensg,
                               colData   = col_data,
                               design    = ~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]   # filter lowly expressed genes
dds <- DESeq(dds)

# Contrast: case vs Control — positive log2FC = up-regulated in case
res    <- results(dds, contrast=c("condition","case","control"))
res_df <- as.data.frame(res)
res_df$ENSEMBL <- rownames(res_df)
res_df <- res_df %>% filter(!is.na(padj))

# Map Ensembl IDs to gene symbols
sym_map <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys    = res_df$ENSEMBL,
                                  columns = c("SYMBOL","GENENAME"),
                                  keytype = "ENSEMBL") %>%
  distinct(ENSEMBL, .keep_all=TRUE)
res_df <- left_join(res_df, sym_map, by="ENSEMBL")

write.table(res_df, "DESeq2_case_vs_Control_RNA.txt",
            quote=FALSE, row.names=FALSE, sep="\t")

n_sig   <- sum(abs(res_df$log2FoldChange) > 1 & res_df$padj < 0.05, na.rm=TRUE)
n_up    <- sum(res_df$log2FoldChange >  1 & res_df$padj < 0.05, na.rm=TRUE)
n_down  <- sum(res_df$log2FoldChange < -1 & res_df$padj < 0.05, na.rm=TRUE)
n_total <- nrow(res_df)

cat("Total genes tested (DESeq2)          :", n_total, "\n")
cat("Significant genes (padj<0.05, |FC|>1):", n_sig,   "\n")
cat("  Up in case                        :", n_up,    "\n")
cat("  Down in case                      :", n_down,  "\n")

# RNA-seq volcano plot
res_df <- res_df %>%
  mutate(regulation = case_when(
    log2FoldChange >  1 & padj < 0.05 ~ "Up in case",
    log2FoldChange < -1 & padj < 0.05 ~ "Down in case",
    TRUE ~ "NS"))

p_volcano <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=regulation)) +
  geom_point(alpha=0.5, size=0.8) +
  scale_color_manual(Treat1ues=c("Up in case"="red","Down in case"="blue","NS"="grey70")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dotted") +
  labs(title="RNA-seq Volcano: case vs Control",
       x="log2FC (case/Control)", y="-log10(padj)") +
  theme_bw()
ggsave("Volcano_RNA_case_vs_Control.pdf", p_volcano, width=7, height=6)

# =============================================================================
# STEP 4: Four-quadrant ATAC-RNA integration
# =============================================================================
cat("\n=== STEP 4: Four-quadrant ATAC vs RNA integration ===\n")

# One represTreat2ative ATAC peak per gene (lowest FDR)
atac_per_gene <- atac %>%
  group_by(SYMBOL) %>%
  slice_min(FDR, n=1, with_ties=FALSE) %>%
  ungroup() %>%
  dplyr::select(SYMBOL, ATAC_Fold=Fold, ATAC_FDR=FDR, distanceToTSS, annotation)

rna_per_gene <- res_df %>%
  filter(!is.na(SYMBOL)) %>%
  arrange(padj) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%   # one row per SYMBOL: keep most significant
  dplyr::select(SYMBOL, RNA_log2FC=log2FoldChange, RNA_padj=padj)

integrated <- inner_join(atac_per_gene, rna_per_gene, by="SYMBOL") %>%
  filter(!is.na(RNA_log2FC), !is.na(ATAC_Fold))

cat("Genes with both ATAC and RNA data:", nrow(integrated), "\n")

integrated <- integrated %>%
  mutate(
    ATAC_sig = ATAC_FDR < 0.01,
    RNA_sig  = RNA_padj < 0.05 & abs(RNA_log2FC) > 1,
    # Negate ATAC_Fold so that the x-axis reads intuitively:
    # negative x = closed in case, positive x = opened in case
    ATAC_accessibility_change = -ATAC_Fold,
    quadrant = case_when(
      ATAC_Fold > 0 & RNA_log2FC > 0 ~ "Q1: ATAC closed & RNA up",
      ATAC_Fold > 0 & RNA_log2FC < 0 ~ "Q2: ATAC closed & RNA down",
      ATAC_Fold < 0 & RNA_log2FC > 0 ~ "Q3: ATAC opened & RNA up",
      ATAC_Fold < 0 & RNA_log2FC < 0 ~ "Q4: ATAC opened & RNA down",
      TRUE ~ "Other"
    )
  )

sig_both_integrated <- integrated %>% filter(ATAC_sig & RNA_sig)

quadrant_counts <- sig_both_integrated %>%
  count(quadrant) %>%
  arrange(desc(n))

cat("\nFour-quadrant counts (genes significant in both assays):\n")
print(quadrant_counts)

n_Q1 <- sum(sig_both_integrated$quadrant == "Q1: ATAC closed & RNA up")
n_Q2 <- sum(sig_both_integrated$quadrant == "Q2: ATAC closed & RNA down")
n_Q3 <- sum(sig_both_integrated$quadrant == "Q3: ATAC opened & RNA up")
n_Q4 <- sum(sig_both_integrated$quadrant == "Q4: ATAC opened & RNA down")
cat("\nQ1 (ATAC closed & RNA up)  :", n_Q1, "genes\n")
cat("Q2 (ATAC closed & RNA down):", n_Q2, "genes\n")
cat("Q3 (ATAC opened & RNA up)  :", n_Q3, "genes\n")
cat("Q4 (ATAC opened & RNA down):", n_Q4, "genes\n")

# Save per-quadrant gene lists
for (q in c("Q1","Q2","Q3","Q4")) {
  q_label <- switch(q,
    Q1 = "Q1: ATAC closed & RNA up",
    Q2 = "Q2: ATAC closed & RNA down",
    Q3 = "Q3: ATAC opened & RNA up",
    Q4 = "Q4: ATAC opened & RNA down"
  )
  q_genes <- sig_both_integrated %>%
    filter(quadrant == q_label) %>%
    dplyr::select(SYMBOL, ATAC_Fold, ATAC_FDR,
                  RNA_log2FC, RNA_padj, distanceToTSS, annotation)
  fname <- paste0("Quadrant_", q, "genes_case2control.txt")
  write.table(q_genes, fname, quote=FALSE, row.names=FALSE, sep="\t")
  cat("Saved:", fname, "\n")
}

write.table(integrated, "ATAC_RNA_integrated_case2control.txt",
            quote=FALSE, row.names=FALSE, sep="\t")

# =============================================================================
# STEP 5: Generate BED file for Q1 closed peaks (for HOMER motif analysis)
# Q1 definition: ATAC closed (Fold > 0) AND RNA up (log2FC > 1, padj < 0.05)
# Input for motif analysis: the actual peak coordinates (not gene coordinates)
# =============================================================================
cat("\n=== STEP 5: Generating Q1 closed peaks BED file for HOMER motif analysis ===\n")

# Get Q1 gene symbols
Q1_symbols <- sig_both_integrated %>%
  filter(quadrant == "Q1: ATAC closed & RNA up") %>%
  pull(SYMBOL)

cat("Q1 genes (ATAC closed & RNA up):", length(Q1_symbols), "\n")

# Retrieve all closed peaks associated with Q1 genes
# Each gene may have multiple closed peaks — include all of them for motif analysis
Q1_closed_peaks <- atac_closed %>%
  filter(SYMBOL %in% Q1_symbols) %>%
  dplyr::select(seqnames, start, end, SYMBOL, FDR, Fold, distanceToTSS)

cat("Total closed peaks associated with Q1 genes:", nrow(Q1_closed_peaks), "\n")
cat("  (Multiple peaks per gene are all included for motif analysis)\n")

# Write BED file (chr, start, end, name, score, strand)
# HOMER requires 0-based BED; DiffBind outputs 1-based — subtract 1 from start
Q1_bed <- data.frame(
  chr    = Q1_closed_peaks$seqnames,
  start  = Q1_closed_peaks$start - 1,   # convert to 0-based
  end    = Q1_closed_peaks$end,
  name   = paste0(Q1_closed_peaks$SYMBOL, "_",
                  Q1_closed_peaks$seqnames, ":",
                  Q1_closed_peaks$start, "-",
                  Q1_closed_peaks$end),
  score  = round(-log10(Q1_closed_peaks$FDR + 1e-300), 2),
  strand = "."
)

write.table(Q1_bed, "Q1_closed_peaks_for_HOMER.bed",
            quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

cat("Saved: Q1_closed_peaks_for_HOMER.bed\n")
cat("  -> Use this file as input for findMotifsGenome.pl\n")
cat("  -> Recommended background: all_closed_peaks.bed (Fold > 0)\n\n")

# Also save all closed peaks BED for use as HOMER background
all_closed_bed <- data.frame(
  chr    = atac_closed$seqnames,
  start  = atac_closed$start - 1,
  end    = atac_closed$end,
  name   = paste0(atac_closed$seqnames, ":",
                  atac_closed$start, "-",
                  atac_closed$end),
  score  = round(-log10(atac_closed$FDR + 1e-300), 2),
  strand = "."
)
write.table(all_closed_bed, "all_closed_peaks_for_HOMER_bg.bed",
            quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
cat("Saved: all_closed_peaks_for_HOMER_bg.bed (for use as HOMER background)\n")

# Scatter plot: ATAC accessibility change vs RNA log2FC
sig_both  <- integrated %>% filter(ATAC_sig, RNA_sig)
top_label <- sig_both %>%
  mutate(score = abs(RNA_log2FC) * abs(ATAC_Fold)) %>%
  slice_max(score, n=30)

p_scatter <- ggplot(integrated,
                    aes(x=ATAC_accessibility_change, y=RNA_log2FC)) +
  # All genes as grey background
  geom_point(color="grey85", size=0.4, alpha=0.6) +
  # Doubly significant genes coloured by quadrant
  geom_point(data=sig_both, aes(color=quadrant), size=1.5, alpha=0.85) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_text_repel(data=top_label, aes(label=SYMBOL),
                  size=2.5, max.overlaps=25, segmTreat2.color="grey50") +
  scale_color_manual(Treat1ues=c(
    "Q1: ATAC closed & RNA up"   = "#E41A1C",
    "Q2: ATAC closed & RNA down" = "#377EB8",
    "Q3: ATAC opened & RNA up"   = "#4DAF4A",
    "Q4: ATAC opened & RNA down" = "#984EA3")) +
  annotate("text", x=-Inf, y= Inf, hjust=-0.1, vjust=1.5,
           label=paste0("Q1 (n=",sum(sig_both$quadrant=="Q1: ATAC closed & RNA up"),")"),
           color="#E41A1C", size=3.5) +
  annotate("text", x= Inf, y= Inf, hjust=1.1,  vjust=1.5,
           label=paste0("Q3 (n=",sum(sig_both$quadrant=="Q3: ATAC opened & RNA up"),")"),
           color="#4DAF4A", size=3.5) +
  annotate("text", x=-Inf, y=-Inf, hjust=-0.1, vjust=-0.5,
           label=paste0("Q2 (n=",sum(sig_both$quadrant=="Q2: ATAC closed & RNA down"),")"),
           color="#377EB8", size=3.5) +
  annotate("text", x= Inf, y=-Inf, hjust=1.1,  vjust=-0.5,
           label=paste0("Q4 (n=",sum(sig_both$quadrant=="Q4: ATAC opened & RNA down"),")"),
           color="#984EA3", size=3.5) +
  labs(
    title = "ATAC-seq vs RNA-seq: case vs Control (ATL43)",
    x     = "ATAC accessibility change\n← Closed in case  |  Opened in case →",
    y     = "RNA log2FC (case/Control)\n↓ Down  |  Up ↑",
    color = "Quadrant"
  ) +
  theme_bw(base_size=12)
ggsave("ATAC_RNA_scatter_case2control.pdf", p_scatter, width=10, height=8)

# =============================================================================
# Save objects needed by PART 2 (clusterProfiler GO enrichmTreat2 steps)
# =============================================================================
cat("\n=== Saving objects for PART 2 (GO enrichmTreat2) ===\n")

save(atac, atac_closed, atac_opened,
     genes_closed, genes_opened,
     res_df, integrated, sig_both_integrated,
     Q1_symbols, Q1_closed_peaks,
     OUTDIR,
     file = file.path(OUTDIR, "ATLL43_integration_part1.RData"))

cat("Saved: ATLL43_integration_part1.RData\n")

# =============================================================================
# Summary
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("PART 1 complete. Summary statistics:\n\n")

cat("--- ATAC-seq ---\n")
cat("  Total diff peaks (FDR<0.01)    :", nrow(atac), "\n")
cat("  Closed peaks (Fold > 0)        :", nrow(atac_closed), "\n")
cat("  Opened peaks (Fold < 0)        :", nrow(atac_opened), "\n")
cat("  Unique genes near closed peaks :", nrow(genes_closed), "\n")
cat("  Unique genes near opened peaks :", nrow(genes_opened), "\n\n")

cat("--- RNA-seq (DESeq2) ---\n")
cat("  Total genes tested             :", n_total, "\n")
cat("  Significant (padj<0.05, |FC|>1):", n_sig,   "\n")
cat("  Up in case                    :", n_up,    "\n")
cat("  Down in case                  :", n_down,  "\n\n")

cat("--- ATAC + RNA Integration ---\n")
cat("  Genes in both datasets         :", nrow(integrated), "\n")
cat("  Q1 (ATAC closed & RNA up)      :", n_Q1, "\n")
cat("  Q2 (ATAC closed & RNA down)    :", n_Q2, "\n")
cat("  Q3 (ATAC opened & RNA up)      :", n_Q3, "\n")
cat("  Q4 (ATAC opened & RNA down)    :", n_Q4, "\n\n")

cat("--- HOMER Motif Analysis Input ---\n")
cat("  Q1 closed peaks (BED)          :", nrow(Q1_bed), "peaks\n")
cat("  Background closed peaks (BED)  :", nrow(all_closed_bed), "peaks\n\n")

cat("Output files:\n")
files <- c(
  "ATAC_closed_genes_case2control.txt      — unique genes near ATAC-closed peaks",
  "ATAC_opened_genes_case2control.txt      — unique genes near ATAC-opened peaks",
  "DESeq2_case_vs_Control_RNA.txt          — RNA-seq DESeq2 full results",
  "Volcano_RNA_case_vs_Control.pdf         — RNA-seq volcano plot",
  "ATAC_RNA_integrated_case2control.txt    — ATAC + RNA integrated table (all genes)",
  "ATAC_RNA_scatter_case2control.pdf       — four-quadrant scatter plot",
  "Quadrant_Q1_genes_case2control.txt      — Q1 gene list (ATAC closed & RNA up)",
  "Quadrant_Q2_genes_case2control.txt      — Q2 gene list (ATAC closed & RNA down)",
  "Quadrant_Q3_genes_case2control.txt      — Q3 gene list (ATAC opened & RNA up)",
  "Quadrant_Q4_genes_case2control.txt      — Q4 gene list (ATAC opened & RNA down)",
  "Q1_closed_peaks_for_HOMER.bed            — Q1 peaks for findMotifsGenome.pl",
  "all_closed_peaks_for_HOMER_bg.bed        — background peaks for HOMER",
  "ATLL43_integration_part1.RData           — objects for PART 2 (GO enrichmTreat2)"
)
cat(paste0("  ", files, "\n"))
cat("\nHOMER command to run next:\n")
cat("  findMotifsGenome.pl Q1_closed_peaks_for_HOMER.bed hg38 ./homer_Q1_closed/ \\\n")
cat("    -size 200 -mask -p 8 \\\n")
cat("    -bg all_closed_peaks_for_HOMER_bg.bed\n")
cat(strrep("=", 60), "\n")
