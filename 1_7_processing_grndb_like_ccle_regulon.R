#NOTE: this script should be run on a server/cluster, as R+RStudio struggle to handle SCENIC
#Input: motif annotation files
#Output: GRNdb-like regulons for use in the downstream analysis

library(tidyverse)
library(SCENIC)
library(RcisTarget)
library(AUCell)

#load dbs
dbs <- list('500bp'='hg19-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather',
            '10kb'='hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather')

#cancer list
cancer_list = c("luad", "coad", "brca", "paad")

for (cancer_label in cancer_list) {
  scenicOptions <- try(initializeScenic(org = "hgnc", dbDir = "./cisTarget_databases", dbs = dbs, nCores = 10))
  
  motifAnnotations_hgnc <- motifAnnotations
  
  scenicOptions <- initializeScenic(org = "hgnc", dbDir = "./cisTarget_databases", dbs = dbs, nCores = 10)
  
  #load gene expression data
  expr_cancer_ccle <- read.csv(paste0("./ccle_regulons_input/", cancer_label, "_expr.txt"), sep = "\t", row.names = "gene") %>%
    as.matrix()
  
  #get genie3 networks
  genesKept <- geneFiltering(expr_cancer_ccle, scenicOptions)
  expr_cancer_ccle_filtered <- expr_cancer_ccle[genesKept, ]
  grn_cancer <- runGenie3(expr_cancer_ccle_filtered, scenicOptions)
  grn_cancer <- grn_cancer %>%
    rename(c("tf" = "TF", "target" = "Target"))
  write.csv(grn_cancer, file = paste0("./",cancer_label, "_genie3_weights.csv"), quote = FALSE, row.names = FALSE) #saving genie3 networks
  
  runCorrelation(expr_cancer_ccle_filtered, scenicOptions)
  
  scenicOptions@settings$dbs <- scenicOptions@settings$dbs
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
  scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
  scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, expr_cancer_ccle_filtered)
  
  #figure out where this goes - check server
  #load NES df for filtering the genie3 grn 
  cancer_regulon_NES_processed <- read.csv("./output/Step2_regulonTargetsInfo.tsv", sep = "\t") %>%
    dplyr::select("TF", "gene")
  
  #load genie3 network #uncomment if needed
  #grn_cancer <- read.csv(paste0(cancer_label, "_genie3_weights.csv"))
  
  #filter genie 3 network based on the edges present in the NES df 
  grn_cancer_filtered <- cancer_regulon_NES_processed %>%
    left_join(grn_cancer, by = c("TF" = "tf", "gene" = "target")) %>%
    na.omit()
  
  colnames(grn_cancer_filtered) <- c("tf", "target", "likelihood") #setting column names
  
  t_expr_cancer_ccle <- expr_cancer_ccle %>%
    t() %>%
    as.data.frame()
  
  #estimate mode of regulation (activation/inhibition) from spearman correlation sign
  for (i in 1:length(grn_cancer_filtered$tf)) { 
    tf <- as.character(grn_cancer_filtered[i, "tf"])
    target <- as.character(grn_cancer_filtered[i, "target"])
    grn_cancer_filtered[i, "mor"] <- cor.test(t_expr_cancer_ccle[, tf], t_expr_cancer_ccle[ ,target], method = "spearman", exact = FALSE)$estimate
  }
  
  #multiply sign of MOR with likelihood to get updated MOR
  grn_cancer_filtered$mor <- sign(grn_cancer_filtered$mor) * grn_cancer_filtered$likelihood 
  grn_cancer_filtered <- grn_cancer_filtered %>% #remove likelihood column as it's not needed by decoupleR
    subset(select = -likelihood)
  
  #save regulons for further use
  write.csv(grn_cancer_filtered, file = paste0("./regulons/grndb_ccle_regulons/", cancer_label, "_grndb_ccle_regulon.csv"), row.names = FALSE)
  print(paste("DONE updating", cancer, "GRNdb regulon"))
  rm(scenicOptions)
  }
