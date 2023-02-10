library(tidyverse)
library(magrittr)
library(stringr)

#choose one: 
regulons <- "aracne" #OR
regulons <- "grndb" #OR
regulons <- "dorothea"

#choose between pancancer and cancer-type analysis
cancer_list = c("aml", "blca", "brca", "coad",
                "gbm", "hnsc", "kirc", "luad",
                "paad", "stad") #OR

cancer_list = c("pancancer") #for pancancer analysis


methods <- c("ULM", "VIPER", "W.Mean", "Consensus", "Expression")


for (cancer_label in cancer_list) {
  cles <- read.csv(paste0("./results/cles/", regulons, "/", cancer_label, "_cles_per_method_per_gene.csv"), row.names = 1)
  p.value <- read.csv(paste0("./results/cles/", regulons, "/", cancer_label, "_p_value_wilcox_per_method_per_gene.csv"), row.names = 1)
  
  #create empty df
  bar_signif_correlations <- data.frame(matrix(, nrow = 3, ncol = length(methods))) 
  
  colnames(bar_signif_correlations) <- methods
  
  #pair each method with CLES values we stratify by
  for (method in methods) {
    bar_signif_correlations[ ,method] <- c("0.7<CLES≤0.8", "0.8<CLES≤0.9", "0.9<CLES≤1")
  }
  
  #pivot to have each pairing in df
  bar_signif_correlations <- bar_signif_correlations %>%
    pivot_longer(cols = everything(), names_to = "Method", values_to = "CLES")
  
  #what cancer type we are using
  bar_signif_correlations$Cancer <- cancer_label
  
  #empty columns to populate
  bar_signif_correlations$signif_drivers <- NA
  bar_signif_correlations$total_drivers <- length(rownames(cles))
  bar_signif_correlations$Percentage_signif <- NA
  
  
  for(threshold in c(0.8, 0.9, 1)){ #thresholds used for CLES
    for (method in methods) {
      lower_threshold <- threshold - 0.1
      bar_signif_correlations[bar_signif_correlations$CLES == paste0(lower_threshold, "<CLES≤", threshold) &
                                bar_signif_correlations$Method == method, "signif_drivers"] <- 
        sum(p.value[ ,method] < 0.05 &
              cles[ ,method] <= threshold &
              cles[ ,method] > lower_threshold, 
            na.rm = TRUE)
    }
  }
  
  #also calculating percentage of significant correlations
  bar_signif_correlations$Percentage_signif <- (bar_signif_correlations$signif_drivers/bar_signif_correlations$total_drivers)*100
  
  write.csv(bar_signif_correlations, 
            file = paste0("./results/cles/", regulons, "/", 
                          cancer_label, "_", regulons, "_number_signif_cles.csv"), row.names = FALSE)
}

