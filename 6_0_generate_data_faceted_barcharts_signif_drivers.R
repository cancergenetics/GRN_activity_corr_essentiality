#This scripts generates the data in order to create Fig 4B,C,D, S3
#Input: correlations between activity and essentiality, no. cell lines per cancer type
#Output: data to be fed for batcharts in script 6_1

library(tidyverse)

#select regulons
regulons <- "aracne" #OR
regulons <- "grndb" #OR
regulons <- "dorothea"

#select cancer list analysis
cancer_list = c("aml", "blca", "brca", "coad",
                "gbm", "hnsc", "kirc", "luad",
                "paad", "stad") #OR
cancer_list = c("pancancer") #only works with dorothea regulons

methods <- c("ULM", "VIPER", "W.Mean", "Consensus", "Expression")

#read files according to type of analysis chosen: pancancer or per cancer type
for (cancer_label in cancer_list) {
  if(cancer_label == "pancancer") {
    if (regulons == "dorothea") { #regulons dorothea, pancancer
      corr_driver_act_exp_matrix <- read.csv("./results/results_matrices_per_combination/dorothea/dorothea_pancancer/1percent_filter_corr_scores_filtered_dorothea_pancancer.csv", row.names = 1)
      
      no_of_cell_lines <- read.csv("./results/results_matrices_per_combination/dorothea/dorothea_pancancer/dorothea_pancancer_dorothea_consensus.csv", row.names = 1) %>%
        ncol()
    } else { #regulons not dorothea, pancancer
           print("You can't pair cancer type-specific regulons with the pancancer analysis")

    }
  } else { 
    if (regulons == "dorothea") {
      corr_driver_act_exp_matrix <- read.csv(paste0("./results/results_matrices_per_combination/", 
                                                    regulons, "/", regulons, "_", cancer_label,
                                                    "/corr_scores_filtered_", regulons, "_", cancer_label, 
                                                    ".csv"), row.names = 1)
      
      no_of_cell_lines <- read.csv(paste0("./results/results_matrices_per_combination/", 
                                          regulons, "/", regulons, "_", cancer_label,
                                          "/", regulons, "_", cancer_label, "_", regulons, "_",
                                          "consensus.csv"), row.names = 1) %>%
        ncol()
    } else { #regulons not dorothea
      corr_driver_act_exp_matrix <- read.csv(paste0("./results/results_matrices_per_combination/", 
                                                    regulons, "/", regulons, "_", cancer_label, "_", cancer_label, 
                                                    "/corr_scores_filtered_", regulons, "_", cancer_label, "_", cancer_label,
                                                    ".csv"), row.names = 1)
      
      no_of_cell_lines <- read.csv(paste0("./results/results_matrices_per_combination/", 
                                          regulons, "/", regulons, "_", cancer_label, "_", cancer_label,
                                          "/", regulons, "_", cancer_label, "_", cancer_label, "_",
                                          "consensus.csv"), row.names = 1) %>%
        ncol()
    }
  }
  
  #1. creating data for absolute values correlations
  #create empty df
  bar_signif_correlations <- data.frame(matrix(, nrow = 4, ncol = length(methods))) 
  colnames(bar_signif_correlations) <- methods
  
  #pair each method with Pearson's R values we stratify by
  for (method in methods) {
    bar_signif_correlations[ ,method] <- c("0.2<|R|≤0.4", "0.4<|R|≤0.6", "0.6<|R|≤0.8", "0.8<|R|≤1")
  }
  
  #pivot to have each pairing in df
  bar_signif_correlations <- bar_signif_correlations %>%
    pivot_longer(cols = everything(), names_to = "Method", values_to = "Pearsons_R")
  
  #what cancer we using
  bar_signif_correlations$Cancer <- cancer_label
  
  #empty columns to populate
  bar_signif_correlations$signif_drivers <- NA
  bar_signif_correlations$total_drivers <- length(corr_driver_act_exp_matrix$hgnc_symbol)
  bar_signif_correlations$Percentage_signif <- NA
  bar_signif_correlations$Number_cell_lines <- no_of_cell_lines
  
  
  for(threshold in c(0.4, 0.6, 0.8, 1)){
    for (method in methods) {
      lower_threshold <- threshold - 0.2
      bar_signif_correlations[bar_signif_correlations$Pearsons_R == paste0(lower_threshold, "<|R|≤", threshold) &
                                bar_signif_correlations$Method == method, "signif_drivers"] <- 
        sum(corr_driver_act_exp_matrix[ ,paste0(method, ".p.value")] < 0.05 &
              abs(corr_driver_act_exp_matrix[ ,paste0(method, ".cor")]) <= threshold &
              abs(corr_driver_act_exp_matrix[ ,paste0(method, ".cor")]) > lower_threshold, 
            na.rm = TRUE)
    }
  }
  
  bar_signif_correlations$Percentage_signif <- (bar_signif_correlations$signif_drivers/bar_signif_correlations$total_drivers)*100
  
  write.csv(bar_signif_correlations, 
            file = paste0("./results/barcharts_data/", 
                          cancer_label, "_", regulons, "_number_signif_corr.csv"), row.names = FALSE)
  
  
  #2. creating data for positive v negative barcharts
  #create empty df
  bar_pos_vs_neg <- data.frame(matrix(, nrow = 2, ncol = length(methods)))
  
  colnames(bar_pos_vs_neg) <- methods
  
  #pair each method with Pearson's R values we stratify by
  for (method in methods) {
    bar_pos_vs_neg[ ,method] <- c("R≤0", "R>0")
  }
  
  #pivot to have each pairing in df
  bar_pos_vs_neg <- bar_pos_vs_neg %>%
    pivot_longer(cols = everything(), names_to = "Method", values_to = "Pearsons_R")
  
  #what cancer we are using
  bar_pos_vs_neg$Cancer <- cancer_label
  
  #empty columns to populate
  bar_pos_vs_neg$count <- NA
  
  #iterate through every method and populate the count for positive and negative correlations, respectively
  for (method in methods) {
    bar_pos_vs_neg[bar_pos_vs_neg$Pearsons_R == "R≤0" &
                     bar_pos_vs_neg$Method == method, "count"] <- 
      sum(corr_driver_act_exp_matrix[ ,paste0(method, ".p.value")] < 0.05 &
            corr_driver_act_exp_matrix[ ,paste0(method, ".cor")] <= 0,
          na.rm = TRUE)
    
    bar_pos_vs_neg[bar_pos_vs_neg$Pearsons_R == "R>0" & 
                     bar_pos_vs_neg$Method == method, "count"] <- 
      sum(corr_driver_act_exp_matrix[ ,paste0(method, ".p.value")] < 0.05 &
            corr_driver_act_exp_matrix[ ,paste0(method, ".cor")] > 0,
          na.rm = TRUE)
  }
  
  
  write.csv(bar_pos_vs_neg, 
            file = paste0("./results/barcharts_data/", 
                          cancer_label, "_", regulons, "_pos_vs_neg_corr.csv"), row.names = FALSE)  
  
}
