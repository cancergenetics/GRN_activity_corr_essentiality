###This script processes grndb regulons to update the MOR so it includes directionality for TF-Target interactions
#Input: unprocessed GRNdb regulons and RNA-Seq matrices for each TCGA cancer type
#Output: processed GRNdb regulons

library(tidyverse)

cancer_list = c("brca","coad", "luad", "paad")

#read expression data after it was processed by the 1_3_preprocess_ccle_regulons_input_expr.R script

brca <- read.table("./ccle_regulons_input/brca_expr.txt", header = TRUE)
luad <- read.table("./ccle_regulons_input/luad_expr.txt", header = TRUE)
coad <- read.table("./ccle_regulons_input/coad_expr.txt", header = TRUE)
paad <- read.table("./ccle_regulons_input/paad_expr.txt", header = TRUE)


for (cancer in cancer_list) {
  #loading cancer type specific expression profile using the string label
  expr_profile <- eval(parse(text = cancer))
  rownames(expr_profile) <- expr_profile$gene #set rownames as genes (annotated with X in each cancer type df)
  expr_profile <- expr_profile %>%
    subset(select = -gene)
  expr_profile_mat <- as.matrix(t(expr_profile)) 
  
  df_expr_profile <- as.data.frame(expr_profile_mat) #get df
  
  #loading in regulons processed by ARACNe-AP from script 1_4_run_aracne_ccle.sh
  grn_matrix <- read.csv(paste0("./aracne_output/aracne_", cancer, "_ccle.txt"), sep = "\t") %>% 
    filter(Regulator %in% colnames(df_expr_profile)) %>% #only keeping TFs present in the CCLE expression matrix
    filter(Target %in% colnames(df_expr_profile)) #only keeping target genes present in the CCLE expression matrix
  
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
  write.csv(grn_matrix, file = paste0("./regulons/aracne_ccle_regulons/", cancer, "_aracne_ccle_regulon.csv"), row.names = FALSE)
  print(paste("DONE updating", cancer, "ARACNe CCLE regulon"))
}   
