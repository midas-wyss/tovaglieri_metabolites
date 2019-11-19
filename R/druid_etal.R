x1 <- results(dds, 
                     contrast = c("Group", "H_E", "M_E"),
                     pAdjustMethod = "fdr",
                     cooksCutoff = FALSE)

x2 <- results(dds, 
              contrast = c("Group", "H_E_AB", "M_E_AB"),
              pAdjustMethod = "fdr",
              cooksCutoff = FALSE)

x3 <- results(dds, 
              contrast = c("Group", "H_E_AB", "H_E"),
              pAdjustMethod = "fdr",
              cooksCutoff = FALSE)


# druid on comparisons
x1 <- data_frame(gene_symbol = xx$gene_id,
                        logFC = as.vector(as.numeric(x1$log2FoldChange)),
                        p_val = as.vector(as.numeric(x1$pvalue)),
                        q_val = as.vector(as.numeric(x1$padj)))

g1 <- x1 %>% dplyr::filter(., q_val < 0.05)
z <- AnnotationDbi::select(x = org.Hs.eg.db, 
                           keys = as.character(g1$gene_symbol), 
                           keytype = "SYMBOL", 
                           columns = "ENTREZID")

# d1 <- run_conddr(conddr_data = "/Volumes/HOME/scripts/r/conddr/data/conddr_data_ndrugs=8422_neffects=41412_saveDate=03082018.RData",
#                      dge_data = cbind(g1$logFC, g1$q_val),
#                      pvalue_thr = 0.05,
#                      log2_fold_thr = 1,
#                      gene_symbols = z$SYMBOL,
#                      gene_entrez = z$ENTREZID,
#                      num_random = 10000,
#                      match_effect = "neg")
# 
# 
# x2 <- data_frame(gene_symbol = xx$gene_id,
#                  logFC = as.vector(as.numeric(x2$log2FoldChange)),
#                  p_val = as.vector(as.numeric(x2$pvalue)),
#                  q_val = as.vector(as.numeric(x2$padj)))
# 
# g2 <- x2 %>% dplyr::filter(., q_val < 0.05)
# z <- AnnotationDbi::select(x = org.Hs.eg.db, 
#                            keys = as.character(g2$gene_symbol), 
#                            keytype = "SYMBOL", 
#                            columns = "ENTREZID")
# 
# d2 <- run_conddr(conddr_data = "/Volumes/HOME/scripts/r/conddr/data/conddr_data_ndrugs=8422_neffects=41412_saveDate=03082018.RData",
#                  dge_data = cbind(g2$logFC, g2$q_val),
#                  pvalue_thr = 0.05,
#                  log2_fold_thr = 1,
#                  gene_symbols = z$SYMBOL,
#                  gene_entrez = z$ENTREZID,
#                  num_random = 10000,
#                  match_effect = "neg")
# 
# 
# x3 <- data_frame(gene_symbol = xx$gene_id,
#                  logFC = as.vector(as.numeric(x3$log2FoldChange)),
#                  p_val = as.vector(as.numeric(x3$pvalue)),
#                  q_val = as.vector(as.numeric(x3$padj)))
# 
# g3 <- x3 %>% dplyr::filter(., q_val < 0.05)
# z <- AnnotationDbi::select(x = org.Hs.eg.db, 
#                            keys = as.character(g3$gene_symbol), 
#                            keytype = "SYMBOL", 
#                            columns = "ENTREZID")
# 
# d3 <- run_conddr(conddr_data = "/Volumes/HOME/scripts/r/conddr/data/conddr_data_ndrugs=8422_neffects=41412_saveDate=03082018.RData",
#                  dge_data = cbind(g3$logFC, g3$q_val),
#                  pvalue_thr = 0.05,
#                  log2_fold_thr = 1,
#                  gene_symbols = z$SYMBOL,
#                  gene_entrez = z$ENTREZID,
#                  num_random = 10000,
#                  match_effect = "neg")

# PLOTS
x1 %>% 
  ggplot() + 
  geom_point(aes(x = logFC, y = -log10(q_val)), color = "gray") + 
  geom_point(data = subset(x1, q_val < 0.01 & logFC > 1), aes(x = logFC, y = -log10(q_val)), color = "red") + 
  geom_point(data = subset(x1, q_val < 0.01 & logFC < -1), aes(x = logFC, y = -log10(q_val)), color = "red") + 
  labs(x = "log2(fold change)", y = "-log10(q)", title = "H_E vs M_E") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 24, color = "black"))

x2 %>% 
  ggplot() + 
  geom_point(aes(x = logFC, y = -log10(q_val)), color = "gray") + 
  geom_point(data = subset(x2, q_val < 0.05 & logFC > 1), aes(x = logFC, y = -log10(q_val)), color = "red") + 
  geom_point(data = subset(x2, q_val < 0.05 & logFC < -1), aes(x = logFC, y = -log10(q_val)), color = "red") + 
  labs(x = "log2(fold change)", y = "-log10(q)", title = "Hmmet + EHEC vs Mmmet + EHEC") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 24, color = "black"),
        axis.title = element_text(size = 24, color = "black"),
        title = element_text(size = 24))
ggsave("results/volcano-heab_meab-txp_09122018.pdf", height = 8, width = 11, dpi = 600)

x3 %>% 
  ggplot() + 
  geom_point(aes(x = logFC, y = -log10(q_val)), color = "gray") + 
  geom_point(data = subset(x3, q_val < 0.05 & logFC > 1), aes(x = logFC, y = -log10(q_val)), color = "red") + 
  geom_point(data = subset(x3, q_val < 0.05 & logFC < -1), aes(x = logFC, y = -log10(q_val)), color = "red") + 
  labs(x = "log2(fold change)", y = "-log10(q)", title = "Hmmet + EHEC vs Hmmet") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 24, color = "black"),
        axis.title = element_text(size = 24, color = "black"),
        title = element_text(size = 24))
ggsave("results/volcano-heab_he-txp_09122018.pdf", height = 8, width = 11, dpi = 600)


# d3 %>%  
#   dplyr::arrange(., desc(score)) %>% 
#   dplyr::filter(., probability_random < 0.001) %>% 
#   ggplot() + 
#   geom_point(aes(x = score, y = fct_reorder(drug_name, score), size = number_matches, color = drug_group)) +
#   labs(x = "score", y = NULL) +
#   theme_minimal() + 
#   theme(axis.text = element_text(size = 12, color = "black"),
#         axis.title = element_text(size = 24, color = "black"))
# ggsave("results/druid-heab_he-txp_09122018.pdf", height = 8, width = 11, dpi = 600)
# 
# 
# d2 %>%  
#   dplyr::arrange(., desc(score)) %>% 
#   dplyr::slice(1:20) %>% 
#   ggplot() + 
#   geom_point(aes(x = score, y = fct_reorder(drug_name, score), size = number_matches, color = drug_group)) +
#   labs(x = "score", y = NULL) +
#   theme_minimal() + 
#   theme(axis.text = element_text(size = 12, color = "black"),
#         axis.title = element_text(size = 24, color = "black"))
# ggsave("results/druid-heab_meab-txp_09122018.pdf", height = 8, width = 11, dpi = 600)