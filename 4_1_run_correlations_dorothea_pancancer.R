###This script calculates correlations for dorothea X pancancer with three different thresholds to filter gene essentiality across the panel of cell lines

library(tidyverse)
library(reshape2)

#select the threshold
threshold <- "1percent" #OR
threshold <- "5percent" #OR
threshold <- "10percent"

regulons_types <- "dorothea"
regulon <- regulons_types

if (threshold == "1percent") {
  no_cell_lines_threshold <- 9 #i.e., only genes that are essential in at least 9 cell lines and non-essential in at least 9 cell lines
} else if(threshold == "5percent") {
  no_cell_lines_threshold <- 48
} else if(threshold == "10percent") {
  no_cell_lines_threshold <- 97
} else {
  print("No threshold selected")
}

cancer_list <- c("pancancer")
cancer_label <- "pancancer"
regulons_list <- c("dorothea")

#create empty matrices to populate at the end of the loops over all regulons and cancer types
ulm_minus_exp_corr <- data.frame(matrix(nrow = length(cancer_list), 
                                        ncol = length(regulons_list)))

rownames(ulm_minus_exp_corr) <- cancer_list
colnames(ulm_minus_exp_corr) <- regulons_list

ulm_corr_pearson <- ulm_minus_exp_corr
ulm_vs_exp_p_values <- ulm_minus_exp_corr

wsum_corr_pearson <- ulm_minus_exp_corr
wsum_minus_exp_corr <- ulm_minus_exp_corr
wsum_vs_exp_p_values <- ulm_minus_exp_corr

viper_corr_pearson <- ulm_minus_exp_corr
viper_minus_exp_corr <- ulm_minus_exp_corr
viper_vs_exp_p_values <- ulm_minus_exp_corr

mlm_corr_pearson <- ulm_minus_exp_corr
mlm_minus_exp_corr <- ulm_minus_exp_corr
mlm_vs_exp_p_values <- ulm_minus_exp_corr

wmean_corr_pearson <- ulm_minus_exp_corr
wmean_minus_exp_corr <- ulm_minus_exp_corr
wmean_vs_exp_p_values <- ulm_minus_exp_corr

consensus_corr_pearson <- ulm_minus_exp_corr
consensus_minus_exp_corr <- ulm_minus_exp_corr
consensus_vs_exp_p_values <- ulm_minus_exp_corr

expr_corr_pearson <- ulm_minus_exp_corr

#separate results per each method
input_dir <- paste0("./results/results_matrices_per_combination/",
                    regulons_types, "/", regulon, "_", cancer_label)

results.W.Mean <- read.csv(paste0(input_dir, "/", regulon, "_", cancer_label, "_wmean.csv"), row.names = 1)
results.W.Sum <- read.csv(paste0(input_dir, "/", regulon, "_", cancer_label, "_wsum.csv"), row.names = 1)
results.ULM <- read.csv(paste0(input_dir, "/", regulon, "_", cancer_label, "_ulm.csv"), row.names = 1)
results.MLM <- read.csv(paste0(input_dir, "/", regulon, "_", cancer_label, "_mlm.csv"), row.names = 1)
results.VIPER <- read.csv(paste0(input_dir, "/", regulon, "_", cancer_label, "_viper.csv"), row.names = 1)
results.Consensus <- read.csv(paste0(input_dir, "/", regulon, "_", cancer_label, "_consensus.csv"), row.names = 1)

t_gene_exp_tmp_subset <- read.csv(paste0(input_dir, "/", regulon, "_", cancer_label, "_expression.csv"), row.names = 1, check.names = FALSE)
gene_effect_tmp_signif_reg_genes <- read.csv(paste0(input_dir, "/", regulon, "_", cancer_label, "_gene_effect.csv"), row.names = 1, check.names = FALSE) %>%
  t() %>%
  as.data.frame()

#same but 9 cell lines
non_essential_genes_per_cancer <- c()
count <- 0
for (gene in colnames(gene_effect_tmp_signif_reg_genes)) {#each gene
  for (cell_line in rownames(gene_effect_tmp_signif_reg_genes)) {#each cell line
    
    if (is.na(gene_effect_tmp_signif_reg_genes[cell_line, gene])) {
      count <- count 
    } else {if (gene_effect_tmp_signif_reg_genes[cell_line, gene] > -0.6)  {
      count <- count + 1 #if a gene is essential in a cell line, count goes up; final count is number of times a gene is essential in the cancer
    }
    }
  }
  if (count > no_cell_lines_threshold) non_essential_genes_per_cancer <- c(non_essential_genes_per_cancer, gene) 
  count <- 0
}

#same but 9 cell lines
essential_genes_per_cancer <- c()
count <- 0
for (gene in colnames(gene_effect_tmp_signif_reg_genes)) {#each gene
  for (cell_line in rownames(gene_effect_tmp_signif_reg_genes)) {#each cell line
    
    if (is.na(gene_effect_tmp_signif_reg_genes[cell_line, gene])) {
      count <- count 
    } else {if (gene_effect_tmp_signif_reg_genes[cell_line, gene] < -0.6)  {
      count <- count + 1 #if a gene is essential in a cell line, count goes up; final count is number of times a gene is essential in the cancer
    }
    }
  }
  if (count >= no_cell_lines_threshold) essential_genes_per_cancer <- c(essential_genes_per_cancer, gene) 
  count <- 0
}

#keep only genes essential in at least one cell line
gene_effect_tmp_signif_reg_genes <- 
  gene_effect_tmp_signif_reg_genes[,(colnames(gene_effect_tmp_signif_reg_genes) %in% essential_genes_per_cancer)]

#keep only genes NON-essential in at least one cell line
gene_effect_tmp_signif_reg_genes <- 
  gene_effect_tmp_signif_reg_genes[,(colnames(gene_effect_tmp_signif_reg_genes) %in% non_essential_genes_per_cancer)]


signif_driver_vector <- colnames(gene_effect_tmp_signif_reg_genes)

corr_per_gene_matrix <- data.frame(matrix(, nrow = length(signif_driver_vector), ncol = 13))
rownames(corr_per_gene_matrix) <- signif_driver_vector
colnames(corr_per_gene_matrix) <- c("hgnc_symbol", 
                                    "ULM.cor", "ULM.p.value",
                                    "W.Sum.cor", "W.Sum.p.value",
                                    "VIPER.cor", "VIPER.p.value",
                                    "W.Mean.cor", "W.Mean.p.value",
                                    "MLM.cor", "MLM.p.value", 
                                    "Consensus.cor", "Consensus.p.value")

#fill matrix with correlations results between driver activity and gene effect
for (gene in signif_driver_vector) {
  W.Mean.cor.tmp <- cor.test(t(gene_effect_tmp_signif_reg_genes)[gene, ], t(results.W.Mean[gene,]), method = "pearson")
  W.Sum.cor.tmp <- cor.test(t(gene_effect_tmp_signif_reg_genes)[gene, ], t(results.W.Sum[gene,]), method = "pearson")
  ULM.cor.tmp <- cor.test(t(gene_effect_tmp_signif_reg_genes)[gene, ], t(results.ULM[gene,]), method = "pearson")
  MLM.cor.tmp <- cor.test(t(gene_effect_tmp_signif_reg_genes)[gene, ], t(results.MLM[gene,]), method = "pearson")
  VIPER.cor.tmp <- cor.test(t(gene_effect_tmp_signif_reg_genes)[gene, ], t(results.VIPER[gene,]), method = "pearson")
  Consensus.cor.tmp <- cor.test(t(gene_effect_tmp_signif_reg_genes)[gene, ], t(results.Consensus[gene,]), method = "pearson")
  
  corr_per_gene_matrix[gene, "hgnc_symbol"] <- gene
  
  corr_per_gene_matrix[gene, "W.Mean.cor"] <- W.Mean.cor.tmp$estimate
  corr_per_gene_matrix[gene, "W.Mean.p.value"] <- W.Mean.cor.tmp$p.value
  
  corr_per_gene_matrix[gene, "W.Sum.cor"] <- W.Sum.cor.tmp$estimate
  corr_per_gene_matrix[gene, "W.Sum.p.value"] <- W.Sum.cor.tmp$p.value
  
  corr_per_gene_matrix[gene, "ULM.cor"] <- ULM.cor.tmp$estimate
  corr_per_gene_matrix[gene, "ULM.p.value"] <- ULM.cor.tmp$p.value
  
  corr_per_gene_matrix[gene, "MLM.cor"] <- MLM.cor.tmp$estimate
  corr_per_gene_matrix[gene, "MLM.p.value"] <- MLM.cor.tmp$p.value
  
  corr_per_gene_matrix[gene, "VIPER.cor"] <- VIPER.cor.tmp$estimate
  corr_per_gene_matrix[gene, "VIPER.p.value"] <- VIPER.cor.tmp$p.value
  
  corr_per_gene_matrix[gene, "Consensus.cor"] <- Consensus.cor.tmp$estimate
  corr_per_gene_matrix[gene, "Consensus.p.value"] <- Consensus.cor.tmp$p.value
  
}
rm(gene)  


#create empty matrix (exp stands for expression)
corr_per_gene_exp_matrix <- data.frame(matrix(, nrow = length(signif_driver_vector), ncol = 3))
rownames(corr_per_gene_exp_matrix) <- signif_driver_vector
colnames(corr_per_gene_exp_matrix) <- c("hgnc_symbol", "Expression.cor", "Expression.p.value")

#fill matrix with correlation results between gene expression and gene effect                                    
for (gene in signif_driver_vector) {
  Expression.cor.tmp <- cor.test(t(gene_effect_tmp_signif_reg_genes)[gene, ], t(t_gene_exp_tmp_subset[gene, ]), method = "pearson")
  corr_per_gene_exp_matrix[gene, "hgnc_symbol"] <- gene
  corr_per_gene_exp_matrix[gene, "Expression.cor"] <- Expression.cor.tmp$estimate
  corr_per_gene_exp_matrix[gene, "Expression.p.value"] <- Expression.cor.tmp$p.value
}
rm(gene)



#make sure same genes are in both activity and expression dfs
corr_per_gene_exp_matrix <- corr_per_gene_exp_matrix %>%
  dplyr::filter(corr_per_gene_exp_matrix$hgnc_symbol %in% rownames(corr_per_gene_matrix))

#merge activity and expression dfs
corr_act_exp_matrix <- merge(corr_per_gene_matrix, corr_per_gene_exp_matrix, by = "hgnc_symbol")

#correlation matrix saved before filtering (needed for barplots)
write.csv(corr_act_exp_matrix, file = paste0(input_dir, "/", threshold, "_filter_corr_scores_filtered_", regulon, "_", cancer_label, ".csv"))

corr_act_exp_matrix[ ,2:15] <- abs(corr_act_exp_matrix[ ,2:15]) #takes absolute values

#mean correlations
viper_corr_pearson[cancer_label, regulon] <- mean(corr_act_exp_matrix$VIPER.cor) 
ulm_corr_pearson[cancer_label, regulon] <- mean(corr_act_exp_matrix$ULM.cor)
mlm_corr_pearson[cancer_label, regulon] <- mean(corr_act_exp_matrix$MLM.cor)
wsum_corr_pearson[cancer_label, regulon] <- mean(corr_act_exp_matrix$W.Sum.cor) 
wmean_corr_pearson[cancer_label, regulon] <- mean(corr_act_exp_matrix$W.Mean.cor)
consensus_corr_pearson[cancer_label, regulon] <- mean(corr_act_exp_matrix$Consensus.cor)
expr_corr_pearson[cancer_label, regulon] <- mean(corr_act_exp_matrix$Expression.cor) 

#calculate differences between activity and expression; calculate p values
#NOTE: The differences between activity and expression are not used later in the analysis, they're just for informative purposes
ulm_minus_exp_corr[cancer_label, regulon] <- mean(corr_act_exp_matrix$ULM.cor) - mean(corr_act_exp_matrix$Expression.cor)
ulm_vs_exp_p_values[cancer_label, regulon] <- wilcox.test(x = corr_act_exp_matrix$ULM.cor, y = corr_act_exp_matrix$Expression.cor, paired = TRUE, alternative = "greater")[["p.value"]]

mlm_minus_exp_corr[cancer_label, regulon] <- mean(corr_act_exp_matrix$MLM.cor) - mean(corr_act_exp_matrix$Expression.cor)
mlm_vs_exp_p_values[cancer_label, regulon] <- wilcox.test(x = corr_act_exp_matrix$MLM.cor, y = corr_act_exp_matrix$Expression.cor, paired = TRUE, alternative = "greater")[["p.value"]]

wsum_minus_exp_corr[cancer_label, regulon] <- mean(corr_act_exp_matrix$W.Sum.cor) - mean(corr_act_exp_matrix$Expression.cor)
wsum_vs_exp_p_values[cancer_label, regulon] <- wilcox.test(x = corr_act_exp_matrix$W.Sum.cor, y = corr_act_exp_matrix$Expression.cor, paired = TRUE, alternative = "greater")[["p.value"]]

viper_minus_exp_corr[cancer_label, regulon] <- mean(corr_act_exp_matrix$VIPER.cor) - mean(corr_act_exp_matrix$Expression.cor)
viper_vs_exp_p_values[cancer_label, regulon] <- wilcox.test(x = corr_act_exp_matrix$VIPER.cor, y = corr_act_exp_matrix$Expression.cor, paired = TRUE, alternative = "greater")[["p.value"]]

wmean_minus_exp_corr[cancer_label, regulon] <- mean(corr_act_exp_matrix$W.Mean.cor) - mean(corr_act_exp_matrix$Expression.cor)
wmean_vs_exp_p_values[cancer_label, regulon] <- wilcox.test(x = corr_act_exp_matrix$W.Mean.cor, y = corr_act_exp_matrix$Expression.cor, paired = TRUE, alternative = "greater")[["p.value"]]

consensus_minus_exp_corr[cancer_label, regulon] <- mean(corr_act_exp_matrix$Consensus.cor) - mean(corr_act_exp_matrix$Expression.cor)
consensus_vs_exp_p_values[cancer_label, regulon] <- wilcox.test(x = corr_act_exp_matrix$Consensus.cor, y = corr_act_exp_matrix$Expression.cor, paired = TRUE, alternative = "greater")[["p.value"]]

#saving matrices containing correlations and p values
dir.create(paste0("./results/correlation_matrices/", regulons_types, "/", threshold, "_filter"))
output_dir <- paste0("./results/correlation_matrices/", regulons_types, "/", threshold, "_filter")

write.csv(ulm_corr_pearson, file = paste0(output_dir, "/ulm_corr_pearson_", regulons_types, ".csv"))
write.csv(ulm_minus_exp_corr, file = paste0(output_dir, "/ulm_minus_exp_corr_", regulons_types, ".csv"))
write.csv(ulm_vs_exp_p_values, file = paste0(output_dir, "/ulm_vs_exp_p_values_", regulons_types, ".csv"))

write.csv(mlm_corr_pearson, file = paste0(output_dir, "/mlm_corr_pearson_", regulons_types, ".csv"))
write.csv(mlm_minus_exp_corr, file = paste0(output_dir, "/mlm_minus_exp_corr_", regulons_types, ".csv"))
write.csv(mlm_vs_exp_p_values, file = paste0(output_dir, "/mlm_vs_exp_p_values_", regulons_types, ".csv"))

write.csv(wsum_corr_pearson, file = paste0(output_dir, "/wsum_corr_pearson_", regulons_types, ".csv"))
write.csv(wsum_minus_exp_corr, file = paste0(output_dir, "/wsum_minus_exp_corr_", regulons_types, ".csv"))
write.csv(wsum_vs_exp_p_values, file = paste0(output_dir, "/wsum_vs_exp_p_values_", regulons_types, ".csv"))

write.csv(viper_corr_pearson, file = paste0(output_dir, "/viper_corr_pearson_", regulons_types, ".csv"))
write.csv(viper_minus_exp_corr, file = paste0(output_dir, "/viper_minus_exp_corr_", regulons_types, ".csv"))
write.csv(viper_vs_exp_p_values, file = paste0(output_dir, "/viper_vs_exp_p_values_", regulons_types, ".csv"))

write.csv(wmean_corr_pearson, file = paste0(output_dir, "/wmean_corr_pearson_", regulons_types, ".csv"))
write.csv(wmean_minus_exp_corr, file = paste0(output_dir, "/wmean_minus_exp_corr_", regulons_types, ".csv"))
write.csv(wmean_vs_exp_p_values, file = paste0(output_dir, "/wmean_vs_exp_p_values_", regulons_types, ".csv"))

write.csv(consensus_corr_pearson, file = paste0(output_dir, "/consensus_corr_pearson_", regulons_types, ".csv"))
write.csv(consensus_minus_exp_corr, file = paste0(output_dir, "/consensus_minus_exp_corr_", regulons_types, ".csv"))
write.csv(consensus_vs_exp_p_values, file = paste0(output_dir, "/consensus_vs_exp_p_values_", regulons_types, ".csv"))

write.csv(expr_corr_pearson, file = paste0(output_dir, "/expr_corr_pearson_", regulons_types, ".csv"))