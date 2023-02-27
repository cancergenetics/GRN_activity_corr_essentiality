#Input: correlations for each cancer type/pancancer
#Input: List of oncogenes downloaded from the cancer gene census
#Output: no files, just prints the p-values of Fisher's exact tests #that verify whether 
#there is an enrichment for oncogenes in the genes with higher activity than expression
#Note: to download the cancer gene census data you will need a COSMIC account
#Go at https://cancer.sanger.ac.uk/census and press "Export: CSV" and save the data as "cgc.csv"

library(tidyverse)

#select regulons:
regulons <- "aracne" #OR
regulons <- "grndb" #OR
regulons <- "dorothea"

cancer_list = c("aml", "gbm", "luad", "brca", "coad",
                "paad", "hnsc", "stad", "blca", "kirc")

methods <- c("ULM", "Viper", "W.Mean", "Consensus", "Expression")


oncogenes <- read.csv("./cgc.csv") %>% 
  dplyr::filter(grepl("oncogene", Role.in.Cancer)) %>% #oncogenes filtering
  dplyr::pull(Gene.Symbol)


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
  rownames(contingency_table) <- c("oncogene", "non_oncogene")
  
  contingency_table["oncogene", "consensus_wins"] <- length(consensus_wins[consensus_wins %in% oncogenes])
  contingency_table["oncogene", "expression_wins"] <- length(consensus_wins[expression_wins %in% oncogenes])
  contingency_table["non_oncogene", "consensus_wins"] <- length(consensus_wins) - length(consensus_wins[consensus_wins %in% oncogenes])
  contingency_table["non_oncogene", "expression_wins"] <- length(expression_wins) - length(consensus_wins[expression_wins %in% oncogenes])
  
  fisher_test_list[[cancer_label]] <- fisher.test(contingency_table)[["p.value"]]
}

#print result of Fisher's exact test
fisher_test_list
