# libraries and functions ----
library(readxl)
library(tidyverse)

# Load data ----
metadata <- read_delim("data/16s_data/2019-02-06_metadata_alessio.txt", delim = "\t")
otu <- read_delim("data/16s_data/2019-02-06_otu_table_alessio.txt", delim = "\t")
tax <- read_delim("data/16s_data/2019-02-06_tax_table_alessio.txt", delim = "\t")

full_df <- merge(x = otu, y = tax, by.x = "otu_id")

# how many otus per sample?
# this will measure observed diversity. this is done including singletons.

x <- as.matrix(otu[, -1])
y <- apply(x, 2, function(s) length(which(s != 0)))

diversity_df <- tibble(sample_id = metadata$name,
                           sample_type = metadata$Type,
                           otu_count = y)

diversity_df %>%
  ggplot() + 
  geom_boxplot(aes(x = sample_type,
                   y = otu_count,
                   fill = sample_type)) +
  scale_fill_manual(breaks = c("Human", "Mouse"),
                    values = c("#CC9966", "#CC0066"),
                    labels = c("Human", "Mouse"),
                    name = "") + 
  labs(x = "Sample Type", y = "Observed Diversity") +
  theme_bw() + 
  theme(axis.title = element_text(size = 36, color = "black"),
        axis.text = element_text(size = 18, color = "black"),
        panel.grid = element_blank(),
        legend.position = "none")
ggsave(filename = paste0("results/", format(Sys.Date(), "%Y-%m-%d"), "_observed-diversity.pdf"),
       height = 8,
       width = 11,
       dpi = 600)


# clean up otus ----
# remove zero count otus
otu2 <- otu
tax2 <- tax

nix1 <- which(rowSums(x) == 0)
if (length(nix1) != 0) {
  otu2 <- otu2[-nix1, ]
  tax2 <- tax2[-nix1, ]
  x <- x[-nix1, ]
}

# summarize to phylum level ----
ap <- unique(tax$kingdom)
ap <- ap[!is.na(ap)]
phy_counts <- matrix(0, nrow = length(ap), ncol = ncol(x))
for (i in seq(1, length(ap))) {
  a1 <- which(tax$kingdom == ap[i])
  a2 <- x[a1, ]
  if (length(a1) > 1) {
    a2 <- colSums(a2)
  }
  phy_counts[i, ] <- a2
}

# summarize to genus level ----
ag <- unique(tax$family)
ag <- ag[!is.na(ag)]
genus_counts <- matrix(0, nrow = length(ag), ncol = ncol(x))
for (i in seq(1, length(ag))) {
  a1 <- which(tax$family == ag[i])
  a2 <- x[a1, ]
  if (length(a1) > 1) {
    a2 <- colSums(a2)
  }
  genus_counts[i, ] <- a2
}

# relative abundances ----
# calculate relative abundances of otus per sample
tots <- colSums(phy_counts)
rel_phy <- matrix(0, nrow = nrow(phy_counts), ncol = ncol(phy_counts))
for (i in seq(1, ncol(phy_counts))) {
  rel_phy[, i] <- phy_counts[, i] / tots[i]
}

tots <- colSums(genus_counts)
rel_genus <- matrix(0, nrow = nrow(genus_counts), ncol = ncol(genus_counts))
for (i in seq(1, ncol(genus_counts))) {
  rel_genus[, i] <- genus_counts[, i] / tots[i]
}


rel_df <- tibble(class = c(rep("phylum", ncol(phy_counts) * nrow(phy_counts)), 
                           rep("genus", ncol(genus_counts) * nrow(genus_counts))),
                 sample_id = c(as.vector(sapply(metadata$name, rep, nrow(phy_counts))),
                               as.vector(sapply(metadata$name, rep, nrow(genus_counts)))),
                 group = c(as.vector(sapply(metadata$Type, rep, nrow(phy_counts))),
                           as.vector(sapply(metadata$Type, rep, nrow(genus_counts)))),
                 taxa_name = c(rep(ap, ncol(phy_counts)),
                               rep(ag, ncol(genus_counts))),
                 abundance = c(as.vector(rel_phy),
                               as.vector(rel_genus)))


rel_df %>%
  dplyr::filter(., class == "phylum") %>%
  ggplot() +
  geom_bar(aes(x = sample_id, y = abundance, fill = taxa_name), color = "black", stat = "identity") + 
  facet_grid(. ~ group, scales = "free") +
  labs(x = "", y = "Relative Abundance") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 36),
        strip.text = element_text(size = 18))
ggsave(filename = paste0("results/", format(Sys.Date(), "%Y-%m-%d"), "_PHYLUM-ABUNDANCE.pdf"),
       height = 8,
       width = 11,
       dpi = 600)

rel_df %>%
  dplyr::filter(., class == "genus") %>%
  ggplot() +
  geom_bar(aes(x = sample_id, y = abundance, fill = taxa_name), color = "black", stat = "identity") + 
  facet_grid(. ~ group, scales = "free") +
  labs(x = "", y = "Relative Abundance") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 36),
        strip.text = element_text(size = 18))
ggsave(filename = paste0("results/", format(Sys.Date(), "%Y-%m-%d"), "_GENUS-ABUNDANCE.pdf"),
       height = 8,
       width = 11,
       dpi = 600)




