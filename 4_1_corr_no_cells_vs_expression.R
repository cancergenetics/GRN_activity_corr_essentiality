#calculates correlation between number of cell lines and activity-vs-essentiality correlation (this is for aracne + consensus)

library(tidyverse)

Consensus <- read.csv("./results/correlation_matrices/aracne/consensus_corr_pearson_aracne.csv", row.names = 1) 
gene_effect <- read.csv("./depmap_input_data/wrangled_input_data/gene_effect_filtered.csv", row.names = 1)

cancer_types <- read.csv("./depmap_input_data/wrangled_input_data/filtered_annotated_sample_info.csv") %>%
  na.omit() %>%
  filter(DepMap_ID %in% rownames(gene_effect)) %>%
  plyr::count("cancer") %>% 
  tibble::column_to_rownames(var = "cancer")

cancer_types$R <- NA
for (cancer_label in rownames(cancer_types)) {
  cancer_types[cancer_label, "R"] <- Consensus[cancer_label, paste0("aracne_", cancer_label)] 
}

#calculate correlation
cor.test(cancer_types$freq, cancer_types$R, method = "pearson")

#get scatterplot
cancer_types %>%
  ggplot(aes(x = freq, y = R)) +
  geom_point() +
  labs(x = "No. cell lines") +
  theme_minimal()
