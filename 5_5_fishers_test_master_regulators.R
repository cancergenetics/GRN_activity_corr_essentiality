#Input: correlations for each cancer type/pancancer
#Input: List of MRs downloaded from Paull et al (2021)
#Output: no files, just prints the p-values of Fisher's exact tests #that verify whether 
        #there is an enrichment for master regulators in the genes with higher activity than expression

library(readxl)
library(tidyverse)

#select regulons:
regulons <- "aracne" #OR
regulons <- "grndb" # OR
regulons <- "dorothea"

cancer_list = c("aml", "gbm", "luad", "brca", "coad",
                "paad", "hnsc", "stad", "blca", "kirc")

methods <- c("ULM", "Viper", "W.Mean", "Consensus", "Expression")

master_regulators <- read_xlsx("/home/cosmin/depmap_analysis/mmc2.xlsx", sheet = 2, col_names = FALSE) %>%
  dplyr::select(-"...2") %>%
  mutate(`...1` = replace(`...1`, `...1` == "laml", "aml"))
   


mrs_list <- list()

for (cancer_label in cancer_list) {
  mrs_list[[cancer_label]] <- master_regulators[master_regulators$...1 == cancer_label, -1] %>%
    pivot_longer(cols = everything(), values_to = "MR") %>%
    drop_na() %>%
    pull("MR") %>%
    unique()
}


fisher_test_list <- list()

for (cancer_label in cancer_list) {
  
  if(regulons == "dorothea"){
    
    corr_driver_act_exp_matrix <- read.csv(paste0("./results/results_matrices_per_combination/", 
                                                  regulons, "/", regulons, "_", cancer_label,
                                                  "/corr_scores_filtered_", regulons, "_", cancer_label, 
                                                  ".csv"), row.names = 1)
    
  } else {
    corr_driver_act_exp_matrix <- read.csv(paste0("./results/results_matrices_per_combination/", 
                                                  regulons, "/", regulons, "_", cancer_label, "_", cancer_label, 
                                                  "/corr_scores_filtered_", regulons, "_", cancer_label, "_", cancer_label,
                                                  ".csv"), row.names = 1)
  }
  
  corr_driver_act_exp_matrix <- corr_driver_act_exp_matrix %>% 
    subset(select = c("hgnc_symbol", "Consensus.cor", "Consensus.p.value", "Expression.cor", "Expression.p.value")) %>%
    mutate(diff = Consensus.cor - Expression.cor)
  
  consensus_wins <- corr_driver_act_exp_matrix %>%
    dplyr::filter(Consensus.p.value < 0.05) %>%
    dplyr::filter(Consensus.cor < Expression.cor) %>%
    dplyr::filter(diff < -0.1) %>% #consensus.cor - expression.cor
    pull("hgnc_symbol")
  
  expression_wins <- corr_driver_act_exp_matrix %>%
    dplyr::filter(Expression.p.value < 0.05) %>%
    dplyr::filter(Consensus.cor > Expression.cor) %>% #negative numbers we're looking at   
    dplyr::filter(diff > 0.1) %>%
    pull("hgnc_symbol")
  
  contingency_table <- data.frame(matrix(nrow = 2, ncol = 2))
  colnames(contingency_table) <- c("consensus_wins", "expression_wins")
  rownames(contingency_table) <- c("MR", "non_MR")
  
  contingency_table["MR", "consensus_wins"] <- length(consensus_wins[consensus_wins %in% mrs_list[[cancer_label]]])
  contingency_table["MR", "expression_wins"] <- length(consensus_wins[expression_wins %in% mrs_list[[cancer_label]]])
  contingency_table["non_MR", "consensus_wins"] <- length(consensus_wins) - length(consensus_wins[consensus_wins %in% mrs_list[[cancer_label]]])
  contingency_table["non_MR", "expression_wins"] <- length(expression_wins) - length(consensus_wins[expression_wins %in% mrs_list[[cancer_label]]])
  
  fisher_test_list[[cancer_label]] <- fisher.test(contingency_table)[["p.value"]]
}

#print result of Fisher's exact test
fisher_test_list
