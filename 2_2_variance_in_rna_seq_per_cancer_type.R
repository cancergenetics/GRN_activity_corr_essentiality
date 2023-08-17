#This script calculates the mean variance in RNA-Seq for each cancer type
#Input: Processed DepMap data on gene effect, gene expression and cell line info
#Output: .csv containing RNA-Seq variance per cancer type

library(tidyverse)
library(reshape2)

filtered_annotated_sample_info <- read.csv("./depmap_input_data/wrangled_input_data/filtered_annotated_sample_info.csv")
gene_effect <- read.csv("./depmap_input_data/wrangled_input_data/gene_effect_filtered.csv", row.names = 1)
t_gene_exp <- read.csv("./depmap_input_data/wrangled_input_data/t_gene_exp_filtered.csv", row.names = 1, check.names = FALSE)

cancer_list = c("aml", "gbm", "luad", "coad", "brca", 
                "paad", "hnsc", "stad", "blca", "kirc")

variance <- NULL

for (cancer_type in cancer_list) {
  cell_lines <- filtered_annotated_sample_info %>%
    dplyr::filter(cancer == cancer_type) %>%
    dplyr::pull("DepMap_ID") %>%
    base::intersect(colnames(t_gene_exp)) %>%
    base::intersect(rownames(gene_effect))

  gene_effect_per_cancer <- gene_effect[cell_lines, ]
  gene_exp_per_cancer <- t_gene_exp[ ,cell_lines]
  
  t_gene_effect_per_cancer <- gene_effect_per_cancer %>% 
    t() %>%
    as.data.frame()
  
  #select genes that are NOT essential in at least one cell line
  non_essential_genes_per_cancer <- t_gene_effect_per_cancer[apply(t_gene_effect_per_cancer, 1, function(x) max(x, na.rm = TRUE) > -0.6), drop = FALSE, ] %>%
    rownames()
  
  #select genes that are essential in at least one cell line
  essential_genes_per_cancer <- t_gene_effect_per_cancer[apply(t_gene_effect_per_cancer, 1, function(x) min(x, na.rm = TRUE) < -0.6), drop = FALSE, ] %>%
    rownames()
    
  #keep only genes essential in at least one cell line and NON-essential in at least one cell line
  t_gene_effect_per_cancer <- 
    t_gene_effect_per_cancer[(rownames(t_gene_effect_per_cancer) %in% 
                                intersect(essential_genes_per_cancer, non_essential_genes_per_cancer)), ]
  
  
  #sometimes essential genes per cancer
  gene_exp_per_cancer <- gene_exp_per_cancer[rownames(t_gene_effect_per_cancer), ]

  variance_per_cancer <- apply(gene_exp_per_cancer, 1, var) %>%
    as.data.frame() %>%
    dplyr::rename(variance_per_gene = "." )
  rownames(variance_per_cancer) <- NULL
  variance_per_cancer$cancer_type  <- cancer_type 
  
  variance <- rbind(variance, variance_per_cancer)  
}

#make cancer type a factor
variance$cancer_type <- factor(variance$cancer_type, levels = cancer_list)

#RNA-Seq variance plot for quick visualization
variance %>%
  ggplot(aes(x = cancer_type, y = variance_per_gene)) +
  geom_boxplot(fill = "#61a5c2", size = 0.6) +
  theme_bw(base_size = 16) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
  theme(legend.position = "none") + 
  coord_cartesian(ylim = c(0, 1.5)) +
  labs(x = "Cancer", y = "RNA variance per gene", title = "RNA expression variance")

#calculate mean variance in each cancer type
mean_variance <- variance %>%
  dplyr::group_by(cancer_type) %>%
  dplyr::summarise_at(vars(variance_per_gene), funs(mean(., na.rm=TRUE)))

#save results for linear modelling
write.csv(mean_variance, file = "./results/rna_seq_mean_variance.csv", row.names = FALSE)
