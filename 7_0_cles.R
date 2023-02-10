library(tidyverse)
library(bmbstats)
library(ggpubr)

###choose 1 of:
regulons <- "aracne" #OR
regulons <- "grndb" #OR
regulons <- "dorothea"

###choose 1 of:
cancer_list <- c("aml", "gbm", "luad", 
                 "coad", "brca", "paad",
                 "hnsc", "stad", "blca", "kirc") #OR

cancer_list <- c("pancancer")

methods <- c("Expression", "Consensus", "ULM", "MLM", "VIPER", "W.Mean", "W.Sum")
okabe_7col <- c("#000000", "#009E73", "#D55E00", "#E69F00", "#56B4E9", "#0072B2", "#CC79A7")

###functions ----

#'*select genes that are NOT essential in at least 3 cell lines*

filter_genes_non_essential_at_least_one_cell_line <- function(df_genes_cell_lines) {
  
  gene_vector <- c()
  count <- 0
  for (gene in rownames(gene_effect_cancer)) {#each gene
    for (cell_line in colnames(gene_effect_cancer)) {#each cell line
      
      if (is.na(gene_effect_cancer[gene, cell_line])) {
        count <- count
      } else {if (gene_effect_cancer[gene, cell_line] > -0.6) {
        count <- count + 1 #if a gene is essential in a cell line, count goes up; final count is number of times a gene is essential in the cancer
      }
      }
    }
    if (count >= 3) gene_vector <- c(gene_vector, gene) 
    count <- 0
  }
  
  df_genes_cell_lines <- 
    df_genes_cell_lines[(rownames(df_genes_cell_lines) %in% gene_vector), ]
}


############################################################
#'*select genes that are essential in at least 3 cell lines*

filter_genes_essential_at_least_one_cell_line <- function(df_genes_cell_lines) {
  gene_vector <- c()
  count <- 0
  for (gene in rownames(gene_effect_cancer)) {#each gene
    for (cell_line in colnames(gene_effect_cancer)) {#each cell line
      
      if (is.na(gene_effect_cancer[gene, cell_line])) {
        count <- count
      } else {if (gene_effect_cancer[gene, cell_line] <= -0.6) {
        count <- count + 1 #if a gene is essential in a cell line, count goes up; final count is number of times a gene is essential in the cancer
      }
      }
    }
    if (count >= 3) gene_vector <- c(gene_vector, gene) 
    count <- 0
  }
  
  df_genes_cell_lines <- 
    df_genes_cell_lines[(rownames(df_genes_cell_lines) %in% gene_vector), ]
}


merge_with_essentiality_df <- function(results_matrix, essentiality_df, method) {
  merged_df <- results_matrix[unique(essentiality_df$Gene), ] %>%
    tibble::rownames_to_column(var = "Gene") %>%   #set gene symbols as a column
    tidyr::pivot_longer(cols = -Gene, names_to = "cell_line", values_to = method) %>%
    dplyr::left_join(essentiality_df, by = c("Gene", "cell_line"))
  return(merged_df)
}


calculate_cles_wilcox_p_value <- function(input, method_chosen) {
  ess <- input %>%
    filter(Method == method_chosen) %>%
    filter(essential == 1) %>%
    pull("score")
  non_ess <- input %>%
    filter(Method == method_chosen) %>%
    filter(essential == 0) %>%
    pull("score")
  cles <- CLES(non_ess, ess)
  p.value <- wilcox.test(non_ess, ess)[["p.value"]]
  return(list(cles, p.value))
}



############################
###start of calculations----
cancer_label <- "aml"
for (cancer_label in cancer_list) {
  
  if (regulons == "dorothea") {
    input_dir <- paste0("./results/results_matrices_per_combination/",
                 regulons, "/", regulons, "_", cancer_label)
    results.W.Mean <- read.csv(paste0(input_dir, "/", regulons, "_", cancer_label, "_", regulons, "_wmean.csv"), row.names = 1)
    results.W.Sum <- read.csv(paste0(input_dir, "/", regulons, "_", cancer_label, "_", regulons, "_wsum.csv"), row.names = 1)
    results.ULM <- read.csv(paste0(input_dir, "/", regulons, "_", cancer_label, "_", regulons, "_ulm.csv"), row.names = 1)
    results.MLM <- read.csv(paste0(input_dir, "/", regulons, "_", cancer_label, "_", regulons, "_mlm.csv"), row.names = 1)
    results.VIPER <- read.csv(paste0(input_dir, "/", regulons, "_", cancer_label, "_", regulons, "_viper.csv"), row.names = 1)
    results.Consensus <- read.csv(paste0(input_dir, "/", regulons, "_", cancer_label, "_", regulons, "_consensus.csv"), row.names = 1)
    t_gene_exp_subset <- read.csv(paste0(input_dir, "/", regulons, "_", cancer_label, "_", regulons, "_expression.csv"), row.names = 1)
    gene_effect_cancer <- read.csv(paste0(input_dir, "/", regulons, "_", cancer_label, "_", regulons, "_gene_effect.csv"), row.names = 1)
    
  } else {
    input_dir <- paste0("./results/results_matrices_per_combination/", 
                 regulons, "/", regulons, "_", cancer_label, "_", cancer_label)
    results.W.Mean <- read.csv(paste0(input_dir, "/", regulons, "_", cancer_label, "_", cancer_label, "_", regulons, "_wmean.csv"), row.names = 1)
    results.W.Sum <- read.csv(paste0(input_dir, "/", regulons, "_", cancer_label, "_", cancer_label, "_", regulons, "_wsum.csv"), row.names = 1)
    results.ULM <- read.csv(paste0(input_dir, "/", regulons, "_", cancer_label, "_", cancer_label, "_", regulons, "_ulm.csv"), row.names = 1)
    results.MLM <- read.csv(paste0(input_dir, "/", regulons, "_", cancer_label, "_", cancer_label, "_", regulons, "_mlm.csv"), row.names = 1)
    results.VIPER <- read.csv(paste0(input_dir, "/", regulons, "_", cancer_label, "_", cancer_label, "_", regulons, "_viper.csv"), row.names = 1)
    results.Consensus <- read.csv(paste0(input_dir, "/", regulons, "_", cancer_label, "_", cancer_label, "_", regulons, "_consensus.csv"), row.names = 1)
    t_gene_exp_subset <- read.csv(paste0(input_dir, "/", regulons, "_", cancer_label, "_", cancer_label, "_", regulons, "_expression.csv"), row.names = 1)
    gene_effect_cancer <- read.csv(paste0(input_dir, "/", regulons, "_", cancer_label, "_", cancer_label, "_", regulons, "_gene_effect.csv"), row.names = 1)
  }
  
  gene_effect_reg_genes_filtered <- gene_effect_cancer %>%
    filter_genes_non_essential_at_least_one_cell_line() %>%   #keep only sometimes essential genes
    filter_genes_essential_at_least_one_cell_line %>% 
    tibble::rownames_to_column(var = "Gene") %>%   #set gene symbols as a column
    tidyr::pivot_longer(cols = -Gene, names_to = "cell_line", values_to = "gene_effect") %>%
    dplyr::mutate(essential = case_when(gene_effect < -0.6 ~ 1,  #1 to mark essential genes
                                        gene_effect > -0.6 ~ 0)) %>% #0 to mark non-essential genes
    dplyr::select(-"gene_effect") #the raw gene effect value is not necessary anymore
  
  #merge each method with the essentiality data
  results.W.Mean_merged_ess <- merge_with_essentiality_df(results.W.Mean, gene_effect_reg_genes_filtered, "W.Mean")
  results.W.Sum_merged_ess <- merge_with_essentiality_df(results.W.Mean, gene_effect_reg_genes_filtered, "W.Sum")
  results.ULM_merged_ess <- merge_with_essentiality_df(results.ULM, gene_effect_reg_genes_filtered, "ULM")
  results.MLM_merged_ess <- merge_with_essentiality_df(results.MLM, gene_effect_reg_genes_filtered, "MLM")
  results.VIPER_merged_ess <- merge_with_essentiality_df(results.VIPER, gene_effect_reg_genes_filtered, "VIPER")
  results.Consensus_merged_ess <- merge_with_essentiality_df(results.Consensus, gene_effect_reg_genes_filtered, "Consensus")
  results.Expr_merged_ess <- merge_with_essentiality_df(t_gene_exp_subset, gene_effect_reg_genes_filtered, "Expression")
  
  #joining all dfs of results into a big one
  cols_to_join <- c("cell_line", "essential", "Gene")
  
  results.all_methods_merged_ess <- left_join(results.Expr_merged_ess, results.W.Mean_merged_ess, by = cols_to_join) %>%
    left_join(results.W.Sum_merged_ess, by = cols_to_join) %>%
    left_join(results.ULM_merged_ess, by = cols_to_join) %>%
    left_join(results.MLM_merged_ess, by = cols_to_join) %>%
    left_join(results.VIPER_merged_ess, by = cols_to_join) %>%
    left_join(results.Consensus_merged_ess, by = cols_to_join)
  
  #creating empty dfs to store results from cles and p value from the wilcoxon test
  cles_per_method_per_gene <- as.data.frame(matrix(nrow = length(unique(results.Expr_merged_ess$Gene)), 
                                                   ncol = length(methods)))
  colnames(cles_per_method_per_gene) <- methods
  rownames(cles_per_method_per_gene) <- unique(results.Expr_merged_ess$Gene)
  
  p_value_per_method_per_gene <- as.data.frame(matrix(nrow = length(unique(results.Expr_merged_ess$Gene)), 
                                                      ncol = length(methods)))
  colnames(p_value_per_method_per_gene) <- methods
  rownames(p_value_per_method_per_gene) <- unique(results.Expr_merged_ess$Gene)
  
  for (gene in unique(results.Expr_merged_ess$Gene)) {
    
    #selecting gene by gene
    per_gene_merged_ess <- results.all_methods_merged_ess %>%
      filter(Gene == gene)
    
    per_gene_merged_ess <- per_gene_merged_ess %>%
      pivot_longer(cols = methods, values_to = "score", names_to = "Method")
    
    for (method in methods) {
      cles_per_method_per_gene[gene, method] <- calculate_cles_wilcox_p_value(per_gene_merged_ess, method)[[1]]
      p_value_per_method_per_gene[gene, method] <- calculate_cles_wilcox_p_value(per_gene_merged_ess, method)[[2]]
    }
  }
  
  write.csv(cles_per_method_per_gene, file = paste0("./results/cles/", regulons, "/", cancer_label, "_cles_per_method_per_gene.csv"))
  write.csv(p_value_per_method_per_gene, file = paste0("./results/cles/", regulons, "/", cancer_label, "_p_value_wilcox_per_method_per_gene.csv"))
}
