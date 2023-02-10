###This script processes grndb regulons

library(tidyverse)

cancer_list = c("aml", "blca", "brca", "coad",
                "gbm", "hnsc", "kirc", "luad",
                "paad", "stad")

#read tcga data after it was processed by the separate_tcga_rna_seq_per_cancer.R script
aml <- read.csv("./tcga_data/aml_tcga_processed.csv") 
gbm <- read.csv("./tcga_data/gbm_tcga_processed.csv")
luad <- read.csv("./tcga_data/luad_tcga_processed.csv")
coad <- read.csv("./tcga_data/coad_tcga_processed.csv")
brca <- read.csv("./tcga_data/brca_tcga_processed.csv")
paad <- read.csv("./tcga_data/paad_tcga_processed.csv")
hnsc <- read.csv("./tcga_data/hnsc_tcga_processed.csv")
stad <- read.csv("./tcga_data/stad_tcga_processed.csv")
blca <- read.csv("./tcga_data/blca_tcga_processed.csv") 
kirc <- read.csv("./tcga_data/kirc_tcga_processed.csv") 


for (cancer in cancer_list) {
  
  #loading cancer type specific expression profile using the string label
  expr_profile <- eval(parse(text = cancer))
  rownames(expr_profile) <- expr_profile$X #set rownames as genes (annotated with X in each cancer type df)
  expr_profile <- expr_profile %>%
    subset(select = -X)
  expr_profile_mat <- as.matrix(t(expr_profile)) 
  
  df_expr_profile <- as.data.frame(expr_profile_mat) #get df
  
  #loading in every regulon downloaded from GRNdb
  grn_matrix <- read.csv(paste0("./regulons/grndb_regulons/", cancer, "_TCGA-regulons.txt"), sep = "\t") %>%
    dplyr::filter(Confidence == "High") %>%
    dplyr::select(c("TF", "gene", "Genie3Weight")) %>%
    na.omit() %>%
    filter(TF %in% colnames(df_expr_profile)) %>% #only keeping TFs present in the TCGA expression matrix
    filter(gene %in% colnames(df_expr_profile)) #only keeping target genes present in the TCGA expression matrix
  
  colnames(grn_matrix) <- c("tf", "target", "likelihood") #setting column names

  #estimate mode of regulation (activation/inhibition) from spearman correlation sign
  for (i in 1:length(grn_matrix$tf)) { 
    tf <- as.character(grn_matrix[i, "tf"])
    target <- as.character(grn_matrix[i, "target"])
    grn_matrix[i, "mor"] <- cor.test(df_expr_profile[, tf], df_expr_profile[ ,target], method = "spearman", exact = FALSE)$estimate
  }
  
  #multiply sign of MOR with likelihood to get updated MOR
  grn_matrix$mor <- sign(grn_matrix$mor) * grn_matrix$likelihood 
  grn_matrix <- grn_matrix %>% #remove likelihood column as it's not needed by decoupleR
    subset(select = -likelihood)
  
  #save regulons for further use
  write.csv(grn_matrix, file = paste0("./regulons/grndb_regulons/", cancer, "_grndb_regulon.csv"), row.names = FALSE)
  print(paste("DONE updating", cancer, "GRNdb regulon"))
}   
