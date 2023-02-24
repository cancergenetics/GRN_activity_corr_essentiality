#This script runs GO enrichment to check whether any biological terms are enriched in genes
#where activity performs better than expression or expression performs better than activity in predicting sensitivity to inhibition
#Input: correlation matrices between activity/expression and sensitivity to inhibition
#Output: two files per each cancer type: terms enriched for genes where expression performs better 
                                    #and terms enriched for genes where activity performs better 

library(WebGestaltR)
library(tidyverse)

#choose one: 
regulons <- "aracne" #OR
regulons <- "grndb" #OR
regulons <- "dorothea"

#choose between pancancer and cancer-type analysis
cancer_list <- c("aml", "gbm", "luad", "coad", "brca", 
                 "paad", "hnsc", "stad", "blca", "kirc")

cancer_list <- c("pancancer")

methods <- c("ULM", "MLM", "W.Sum", "W.Mean", "Viper", "Consensus")


for (cancer_type in cancer_list) {
  
  if(cancer_label == "pancancer") {
    if (regulons == "dorothea") { #regulons dorothea, pancancer
      corr_driver_act_exp_matrix <- read.csv("./results/results_matrices_per_combination/dorothea/dorothea_pancancer/1percent_filter_corr_scores_filtered_dorothea_pancancer.csv", row.names = 1)  %>%
        subset(select = c("hgnc_symbol", "Consensus.cor", "Expression.cor"))
    
    } else { #regulons not dorothea, pancancer
      print("You can't pair cancer type-specific regulons with the pancancer analysis")
      
    }
  } else { 
    if (regulons == "dorothea") {
      corr_driver_act_exp_matrix <- read.csv(paste0("./results/results_matrices_per_combination/", 
                                                    regulons, "/", regulons, "_", cancer_label,
                                                    "/corr_scores_filtered_", regulons, "_", cancer_label, 
                                                    ".csv"), row.names = 1) %>%
        subset(select = c("hgnc_symbol", "Consensus.cor", "Expression.cor"))
      

    } else { #regulons not dorothea
      corr_driver_act_exp_matrix <- read.csv(paste0("./results/results_matrices_per_combination/", 
                                                    regulons, "/", regulons, "_", cancer_label, "_", cancer_label, 
                                                    "/corr_scores_filtered_", regulons, "_", cancer_label, "_", cancer_label,
                                                    ".csv"), row.names = 1) %>%
        subset(select = c("hgnc_symbol", "Consensus.cor", "Expression.cor"))
      
    }
  }

  corr_driver_act_exp_matrix$diff_consensus_expression <- abs(corr_driver_act_exp_matrix$Consensus.cor) -
    abs(corr_driver_act_exp_matrix$Expression.cor)
  
  corr_driver_act_exp_matrix <- corr_driver_act_exp_matrix[order(corr_driver_act_exp_matrix$diff_consensus_expression), ]
  
  background_list_cancer <- corr_driver_act_exp_matrix %>%
    dplyr::pull("hgnc_symbol")
  
  consensus_wins_cancer <- corr_driver_act_exp_matrix %>%
    dplyr::filter(diff_consensus_expression > 0) %>%
    dplyr::pull("hgnc_symbol")
  
  expression_wins_cancer <- corr_driver_act_exp_matrix %>%
    dplyr::filter(diff_consensus_expression < 0) %>%
    dplyr::pull("hgnc_symbol")
  
  
  go_expression_wins <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                       "geneontology_Cellular_Component_noRedundant",
                                                       "geneontology_Molecular_Function_noRedundant"),
                                    interestGeneType = "genesymbol",
                                    referenceGeneType = "genesymbol",
                                    interestGene = expression_wins_cancer,
                                    referenceGene = background_list_cancer,
                                    sigMethod = "fdr", fdrThr = 0.1)
  
  go_consensus_wins <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                      "geneontology_Cellular_Component_noRedundant",
                                                      "geneontology_Molecular_Function_noRedundant"),
                                   interestGeneType = "genesymbol",
                                   referenceGeneType = "genesymbol",
                                   interestGene = consensus_wins_cancer,
                                   referenceGene = background_list_cancer,
                                   sigMethod = "fdr", fdrThr = 0.1)
  
  #for our analysis these files will come out empty because there are no enriched terms
  #change the fdrThr parameter to a different threshold to find the "highest" enriched terms. 
     #However, these won't be significant
  write.csv(go_expression_wins, file = paste0("./results/go_enrichment/", regulons, "_", cancer_type, "_expression_wins.csv"), row.names = FALSE, quote = FALSE)
  write.csv(go_consensus_wins, file = paste0("./results/go_enrichment/", regulons, "_", cancer_type, "_consensus_wins.csv"), row.names = FALSE, quote = FALSE)
  
}