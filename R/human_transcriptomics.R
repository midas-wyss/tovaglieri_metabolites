# bioconductor libraries
library(GO.db)
library(reactome.db)
library(DESeq2)
# library(limma)
# library(org.EcK12.eg.db)

# other libraries
library(tidyverse)
library(tidytext)
library(Matrix)
library(gplots)
library(RColorBrewer)
library(readxl)
library(DESeq2)
library(org.Hs.eg.db)
library(DRUID)


# some functions
source("R/enrichment_functions.R")

# rna-seq data 
load(file = "data/human-centric_dual_rnaseq_05162018.RData") # <-- human centric data

samples <- hs_data[[3]]
human_counts <- as.matrix(hs_data[[1]])
human_genes <- hs_data[[2]]

# filter out genes for which we have no gene name
nix <- which(is.na(human_genes$symbol))
if (length(nix) != 0) {
  human_genes <- human_genes[-nix, ]
  human_counts <- human_counts[-nix, ]
}


# remove genes with zero counts
nix <- which(rowSums(human_counts) == 0)

if(length(nix) != 0)
{
  human_counts <- human_counts[-nix, ]
  human_genes <- human_genes[-nix, ]
}


# summarize gene ids
xx <- gene_summarizing(expression_data = human_counts, 
                       gene_id = human_genes$symbol)


# DESeq2 ----
# E. coli
dds <- DESeqDataSetFromMatrix(countData = xx$summarized_expression[,c(1:8, 12:19)],
                              colData = samples[c(1:8, 12:19),],
                              design = ~ Group)

dds <- DESeq(object = dds)

# ec_e_ab_res <- results(dds,contrast=c("Group","H_E_AB","M_E_AB"),alpha=0.05,lfcThreshold = 1)
hc_ab_res <- results(dds, 
                     contrast = c("Group", "H_E_AB", "M_E_AB"),
                     pAdjustMethod = "fdr",
                     cooksCutoff = FALSE)

hc_ab_res <- data_frame(gene_symbol = xx$gene_id,
                        logFC = as.vector(as.numeric(hc_ab_res$log2FoldChange)),
                        p_val = as.vector(as.numeric(hc_ab_res$pvalue)),
                        q_val = as.vector(as.numeric(hc_ab_res$padj)))

diff_genes <- hc_ab_res %>% dplyr::filter(., q_val < 0.05)

# get entrez ids for genes
z <- AnnotationDbi::select(x = org.Hs.eg.db, 
                           keys = as.character(diff_genes$gene_symbol), 
                           keytype = "SYMBOL", 
                           columns = "ENTREZID")

hs_inf <- DRUID::concoct(dge_matrix = cbind(diff_genes$logFC, diff_genes$q_val), 
                         tfidf_matrix = cmap_druid$tfidf, 
                         tfidf_crossproduct = cmap_druid$cpm, 
                         num_random = 10000, 
                         druid_direction = "neg", 
                         fold_thr = 1, 
                         pvalue_thr = 0.05, 
                         entrez = z$ENTREZID)

hs_inf <- hs_inf %>% 
  tibble::add_column(., drug_name = cmap_druid$drugs$name, .before = 1) %>%
  tibble::add_column(., concentration = cmap_druid$drugs$concentration, .before = 2) %>%
  tibble::add_column(., cell_line = cmap_druid$drugs$cell_line, .before = 3)