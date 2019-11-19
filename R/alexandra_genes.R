alexandra_genes <- c("KCNQ1", "KCNQ2", "KCNQ3", "KCNQ4", "KCNQ5", 
                     "CFTR", "PTGER4", "SLC12A2")


rel_abun <- apply(xx$summarized_expression, 2, function(x) x/sum(x))

data_alexandra <- tibble(gene = c(rep(as.character(xx$gene_id[which(xx$gene_id %in% alexandra_genes)]), 4), 
                                  rep(as.character(xx$gene_id[which(xx$gene_id %in% alexandra_genes)]), 4)),
                         group = c(rep("count_data", 4 * length(which(xx$gene_id %in% alexandra_genes))), 
                                   rep("rel_abun", 4 * length(which(xx$gene_id %in% alexandra_genes)))),
                         replicate = c(rep(1, length(which(xx$gene_id %in% alexandra_genes))),
                                       rep(2, length(which(xx$gene_id %in% alexandra_genes))),
                                       rep(3, length(which(xx$gene_id %in% alexandra_genes))),
                                       rep(4, length(which(xx$gene_id %in% alexandra_genes))),
                                       rep(1, length(which(xx$gene_id %in% alexandra_genes))),
                                       rep(2, length(which(xx$gene_id %in% alexandra_genes))),
                                       rep(3, length(which(xx$gene_id %in% alexandra_genes))),
                                       rep(4, length(which(xx$gene_id %in% alexandra_genes)))),
                         data = c(as.vector(xx$summarized_expression[which(xx$gene_id %in% alexandra_genes), 1:4]), 
                                  as.vector(rel_abun[which(xx$gene_id %in% alexandra_genes), 1:4])))

data_alexandra %>% 
  dplyr::filter(., group == "count_data") %>%
  ggplot() + 
  geom_boxplot(aes(x = gene, y = log10(data)), fill = "tan") +
  labs(x = "Gene", y = "log10(summarized counts)") +
  # facet_grid(. ~ group) +
  theme_bw() + 
  theme(axis.text.y = element_text(size = 24, color = "black"),
        axis.text.x = element_text(size = 24, color = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 36, color = "black"),
        axis.ticks = element_line(color = "black"))
ggsave(file = "genes_alexandra.pdf",
       height = 8,
       width = 11,
       dpi = 600)



