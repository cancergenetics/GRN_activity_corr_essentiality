###this scripts gets VIPER, W.Mean, W.Sum, ULM, MLM and Consensus activity values in ten cancer types for each gene
###works with aracne, grndb, and dorothea regulons
###outputs are in separate files in separate directories for each regulon+cancer combination

library(tidyverse)
library(decoupleR)
library(binilib)
library(gplots)
library(org.Hs.eg.db)
library(ggpubr)
library(reshape2)

# loading data ----
filtered_annotated_sample_info <- read.csv("./depmap_input_data/wrangled_input_data/filtered_annotated_sample_info.csv")
gene_effect <- read.csv("./depmap_input_data/wrangled_input_data/gene_effect_filtered.csv", row.names = 1)
t_gene_exp <- read.csv("./depmap_input_data/wrangled_input_data/t_gene_exp_filtered.csv", row.names = 1, check.names = FALSE)

# First select which regulons to use: ----
regulons_types <- "aracne" #OR
regulons_types <- "grndb" #OR
regulons_types <- "dorothea"

cancer_list = c("aml", "gbm", "luad", "coad", "brca", 
                "paad", "hnsc", "stad", "blca", "kirc")

# Functions ----
#separates results from decoupler for each method - !specify method as a string ie. "viper"
sep_decoupler_results <- function(results, method) {
  results.method <- results %>%
    dplyr::filter(statistic == method) %>%
    dplyr::select("source", "condition", "score") %>%
    pivot_wider(names_from = condition, values_from = score) %>%
    as.data.frame()
  rownames(results.method) <- results.method$source
  results.method$source <- NULL
  return(results.method)
}    


# Load regulons ----
if (regulons_types == "aracne") {
  
  #read aracne regulons
  aracne_aml <- read.csv("./regulons/aracne_regulons/aml_aracne_regulon.csv")
  aracne_gbm <- read.csv("./regulons/aracne_regulons/gbm_aracne_regulon.csv")
  aracne_luad <- read.csv("./regulons/aracne_regulons/luad_aracne_regulon.csv")
  aracne_coad <- read.csv("./regulons/aracne_regulons/coad_aracne_regulon.csv")
  aracne_brca <- read.csv("./regulons/aracne_regulons/brca_aracne_regulon.csv")
  aracne_paad <- read.csv("./regulons/aracne_regulons/paad_aracne_regulon.csv")
  aracne_hnsc <- read.csv("./regulons/aracne_regulons/hnsc_aracne_regulon.csv")
  aracne_stad <- read.csv("./regulons/aracne_regulons/stad_aracne_regulon.csv")
  aracne_blca <- read.csv("./regulons/aracne_regulons/blca_aracne_regulon.csv")
  aracne_kirc <- read.csv("./regulons/aracne_regulons/kirc_aracne_regulon.csv")

  regulons_list = c("aracne_aml", "aracne_gbm", "aracne_luad", "aracne_coad", 
                    "aracne_brca", "aracne_paad", "aracne_hnsc", 
                    "aracne_stad", "aracne_blca", "aracne_kirc")
  
} else if(regulons_types == "grndb") {
  
  #read grndb regulons
  grndb_aml <- read.csv("./regulons/grndb_regulons/aml_grndb_regulon.csv")
  grndb_gbm <- read.csv("./regulons/grndb_regulons/gbm_grndb_regulon.csv")
  grndb_luad <- read.csv("./regulons/grndb_regulons/luad_grndb_regulon.csv")          
  grndb_coad <- read.csv("./regulons/grndb_regulons/coad_grndb_regulon.csv")
  grndb_brca <- read.csv("./regulons/grndb_regulons/brca_grndb_regulon.csv")        
  grndb_paad <- read.csv("./regulons/grndb_regulons/paad_grndb_regulon.csv")
  grndb_hnsc <- read.csv("./regulons/grndb_regulons/hnsc_grndb_regulon.csv")
  grndb_stad <- read.csv("./regulons/grndb_regulons/stad_grndb_regulon.csv")
  grndb_blca <- read.csv("./regulons/grndb_regulons/blca_grndb_regulon.csv")
  grndb_kirc <- read.csv("./regulons/grndb_regulons/kirc_grndb_regulon.csv")
  
  
  regulons_list = c("grndb_aml", "grndb_gbm", "grndb_luad", 
                    "grndb_coad", "grndb_brca", 
                    "grndb_paad", "grndb_hnsc", "grndb_stad", 
                    "grndb_blca", "grndb_kirc")

} else if(regulons_types == "dorothea") {
  
  #load dorothea regulons
  data("dorothea_hs", package = "dorothea")
  dorothea = dorothea_hs
  dorothea <- dorothea_hs %>%
    as.data.frame() %>%
    filter(dorothea$tf %in% rownames(t_gene_exp)) %>%
    filter(confidence %in% c("A", "B", "C")) %>%
    subset(select = -confidence)
  
  regulons_list <- c("dorothea")
  
} else {
  print("No regulons selected")
}


# Calculating activity ----
for (cancer_label in cancer_list) { #loop1 for each cancer
  for (regulon in regulons_list){ #loop2 for each regulon
    
    #depmap_id_labels, gene_effect_tmp, t_gene_exp_tmp and regulon.df get updated accordingly with each iteration of the for loop
    
    depmap_id_labels <- filtered_annotated_sample_info %>% #selecting depmap IDs only for specific cancer
      dplyr::filter(cancer == cancer_label) %>%
      dplyr::pull("DepMap_ID")
    
    #cell lines scores for each cancer
    gene_effect_tmp <- gene_effect[rownames(gene_effect) %in% depmap_id_labels, ]
    
    #expression data for each cancer
    t_gene_exp_tmp <- t_gene_exp %>%
      subset(select = rownames(gene_effect_tmp))
    
    t_gene_exp_tmp[t_gene_exp_tmp == 0] <- NA  #replace 0's with NAs
    t_gene_exp_tmp <- t_gene_exp_tmp[which(rowMeans(!is.na(t_gene_exp_tmp)) > 0.8), ] #delete genes with more than 80% of the values NAs
    t_gene_exp_tmp[is.na(t_gene_exp_tmp)] = 0 #replace back NAs with 0s
    
    t_gene_exp_tmp_mat <- as.matrix(t_gene_exp_tmp) #transform into matrix
    
    regulon.df <- as.data.frame(eval(parse(text = regulon))) #get regulon as dataframe
    
    regulon.tmp.target.df <- regulon.df %>%
      dplyr::filter(regulon.df$tf %in% row.names(t_gene_exp_tmp)) %>%
      tidyr::drop_na() #filter to keep only tfs present in expression matrix, otherwise decoupleR throws an error
    
    print(paste("start decoupling", regulon, "+", cancer_label)) #shows progress on the for loop
    
    #decoupleR main chunk----
    results.merged <- decouple(
      mat = t_gene_exp_tmp_mat,
      network = regulon.tmp.target.df,
      .source = "tf",
      .target = "target",
      statistics = c("wmean", "wsum", "ulm", "mlm", "viper"),
      consensus_score = TRUE,
      include_time = FALSE,
      show_toy_call = FALSE,
      args = list(
        wmean = list(.mor = "mor"),
        wsum = list(.mor = "mor"),
        ulm = list(.mor = "mor"),
        mlm = list(.mor = "mor"),
        viper = list(.mor = "mor", pleiotropy = FALSE, verbose = FALSE)
      )
    )
    
    print(paste("Decoupled", regulon,"+", cancer_label)) #shows progress on the for loop
    
    #separate results for each method
    results.VIPER <- sep_decoupler_results(results.merged, "viper")
    results.ULM <- sep_decoupler_results(results.merged, "ulm")
    results.MLM <- sep_decoupler_results(results.merged, "mlm")
    results.W.Sum <- sep_decoupler_results(results.merged, "wsum")
    results.W.Mean <- sep_decoupler_results(results.merged, "wmean")
    results.Consensus <- sep_decoupler_results(results.merged, "consensus")
    
    #makes sure only regulatory genes in all methods are used downstream
    common_genes <- Reduce(intersect,  
                             list(rownames(results.VIPER),rownames(results.ULM), 
                                  rownames(results.W.Sum), rownames(results.W.Mean),
                                  rownames(results.MLM), rownames(results.Consensus)))
    
    gene_effect_tmp_signif_reg_genes <- gene_effect_tmp %>%
      subset(select = common_genes) 
    
    t_gene_exp_tmp_subset <- t_gene_exp_tmp %>%  #makes sure only regulatory genes in all methods are used downstream
      dplyr::filter(rownames(t_gene_exp_tmp) %in% common_genes)
    
    #transpose gene effect dataframe
    gene_effect_tmp_signif_reg_genes <- as.data.frame(t(as.matrix(gene_effect_tmp_signif_reg_genes)))
    
    #create directories for output
    dir.create(paste0("./results/results_matrices_per_combination/", regulons_types, "/", regulon, "_", cancer_label))
    output_dir <- paste0("./results/results_matrices_per_combination/", 
                                regulons_types, "/", regulon, "_", cancer_label)
    
    #save the results
    write.csv(results.VIPER, file = paste0(output_dir, "/", regulon, "_", cancer_label, "_viper.csv"))
    write.csv(results.ULM, file = paste0(output_dir, "/", regulon, "_", cancer_label, "_ulm.csv"))
    write.csv(results.MLM, file = paste0(output_dir, "/", regulon, "_", cancer_label, "_mlm.csv"))
    write.csv(results.W.Sum, file = paste0(output_dir, "/", regulon, "_", cancer_label, "_wsum.csv"))
    write.csv(results.W.Mean, file = paste0(output_dir, "/", regulon, "_", cancer_label, "_wmean.csv"))
    write.csv(results.Consensus, file = paste0(output_dir, "/", regulon, "_", cancer_label, "_consensus.csv"))
    write.csv(t_gene_exp_tmp_subset, file = paste0(output_dir, "/", regulon, "_", cancer_label, "_expression.csv"))
    write.csv(gene_effect_tmp_signif_reg_genes, file = paste0(output_dir, "/", regulon, "_", cancer_label, "_gene_effect.csv"))
    
    print(paste(regulon,"+", cancer_label, "DONE")) #shows progress on the for loop
  }
}
