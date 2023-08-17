##This script pulls individual correlations for genes with at least one R value (expression or activity) < -0.6
#Input: results_matrices_per_combination contents
#Output: Individual scatter plots with line of best fit and Pearsons's R and p-value 
#Output: Supplementary Fig S6B, C

library(stringr)
library(tidyverse)
library(ggpubr)
library(ggrepel)

#select regulons: 
regulons <- "aracne" #OR
regulons <- "grndb" #OR
regulons <- "dorothea"

cancer_list = c("aml", "blca", "brca", "coad",
                "gbm", "hnsc", "kirc", "luad",
                "paad", "stad")

regulons_types <- regulons

methods <- c("ULM", "MLM", "Viper", "W.Mean", "W.Sum", "Consensus", "Expression")


#regulons dictionaries
regulon_cancer_grndb_dict <- c("aml" = "grndb_aml", "gbm" = "grndb_gbm", "luad" = "grndb_luad", 
                               "coad" = "grndb_coad", "brca" = "grndb_brca", "paad" = "grndb_paad", 
                               "hnsc" = "grndb_hnsc", "stad" = "grndb_stad", "blca" = "grndb_blca", 
                               "kirc" = "grndb_kirc")

regulon_cancer_aracne_dict <- c("aml" = "aracne_aml", "gbm" = "aracne_gbm", "luad" = "aracne_luad", 
                                "coad" = "aracne_coad", "brca" = "aracne_brca", "paad" = "aracne_paad", 
                                "hnsc" = "aracne_hnsc", "stad" = "aracne_stad", "blca" = "aracne_blca", 
                                "kirc" = "aracne_kirc")

regulon_cancer_dorothea_dict <- c("aml" = "dorothea", "gbm" = "dorothea", "luad" = "dorothea", 
                                  "coad" = "dorothea", "brca" = "dorothea", "paad" = "dorothea", 
                                  "hnsc" = "dorothea", "stad" = "dorothea", "blca" = "dorothea", 
                                  "kirc" = "dorothea")

signif_genes_method <- NULL

for (cancer_label in cancer_list) {

  if (regulons == "dorothea") {
    corr_driver_act_exp_matrix <- read.csv(paste0("./results/results_matrices_per_combination/", 
                                                  regulons, "/", regulons, "_", cancer_label,
                                                  "/corr_scores_filtered_", regulons, "_", cancer_label, 
                                                  ".csv"), row.names = 1)

      ncol()
  } else { #regulons not dorothea
    corr_driver_act_exp_matrix <- read.csv(paste0("./results/results_matrices_per_combination/", 
                                                  regulons, "/", regulons, "_", cancer_label, "_", cancer_label, 
                                                  "/corr_scores_filtered_", regulons, "_", cancer_label, "_", cancer_label,
                                                  ".csv"), row.names = 1)
  }
  
#only selects genes where ar least one of the methods finds a correlation < -0.6
  signif_genes_method_tmp <- corr_driver_act_exp_matrix %>%
    dplyr::filter(ULM.cor < -0.6 | MLM.cor < -0.6 | VIPER.cor < -0.6 | 
             W.Mean.cor < -0.6 | W.Sum.cor < -0.6 | Consensus.cor < -0.6 | Expression.cor < -0.6) %>%
    subset(select = c("hgnc_symbol", "ULM.cor", "MLM.cor", "VIPER.cor", 
                      "W.Mean.cor", "W.Sum.cor", "Consensus.cor", "Expression.cor"))
  
  if(nrow(signif_genes_method_tmp) > 0) {
  signif_genes_method_tmp$cancer <- cancer_label 
  signif_genes_method <- rbind(signif_genes_method, signif_genes_method_tmp)
  }
}

for (cancer_label in cancer_list)  {
  
  #selects which dictionary to use based on regulons selected
  if (regulons == "grndb") {
    dict <- regulon_cancer_grndb_dict
  } else if (regulons == "aracne") {
    dict <- regulon_cancer_aracne_dict
  } else if(regulons == "dorothea"){
    dict <- regulon_cancer_dorothea_dict
  } else {
    print("No regulons given")
  }
  
  input_dir <- paste0("/home/cosmin/depmap_analysis/results/results_matrices_per_combination/",
                      regulons, "/", dict[cancer_label], "_", cancer_label)
  
  W.Mean <- read.csv(paste0(input_dir, "/", dict[cancer_label], "_", cancer_label, "_", regulons, "_wmean.csv"), row.names = 1, check.names = FALSE)
  W.Sum <- read.csv(paste0(input_dir, "/", dict[cancer_label], "_", cancer_label, "_", regulons, "_wsum.csv"), row.names = 1, check.names = FALSE)
  ULM <- read.csv(paste0(input_dir, "/", dict[cancer_label], "_", cancer_label, "_", regulons, "_ulm.csv"), row.names = 1, check.names = FALSE)
  MLM <- read.csv(paste0(input_dir, "/", dict[cancer_label], "_", cancer_label, "_", regulons, "_mlm.csv"), row.names = 1, check.names = FALSE)
  Viper <- read.csv(paste0(input_dir, "/", dict[cancer_label], "_", cancer_label, "_", regulons, "_viper.csv"), row.names = 1, check.names = FALSE)
  Consensus <- read.csv(paste0(input_dir, "/", dict[cancer_label], "_", cancer_label, "_", regulons, "_consensus.csv"), row.names = 1, check.names = FALSE)
  Expression <- read.csv(paste0(input_dir, "/", dict[cancer_label], "_", cancer_label, "_", regulons, "_expression.csv"), row.names = 1, check.names = FALSE)
  gene_effect_cancer <- read.csv(paste0(input_dir, "/", dict[cancer_label], "_", cancer_label, "_", regulons, "_gene_effect.csv"), row.names = 1, check.names = FALSE)
  
  corr_cancer <- signif_genes_method %>%
    dplyr::filter(cancer == cancer_label)

  for (gene in corr_cancer$hgnc_symbol) {
    df_per_gene <- data.frame(matrix(, nrow = ncol(gene_effect_cancer), ncol = 7))
    colnames(df_per_gene) <- c("ULM", "MLM", "Viper", "W.Mean", "W.Sum", "Consensus", "Expression")
    rownames(df_per_gene) <- colnames(gene_effect_cancer)
    
    for(cell_line in colnames(gene_effect_cancer)) {
      df_per_gene[cell_line, "W.Mean"] <- W.Mean[gene, cell_line]
      df_per_gene[cell_line, "W.Sum"] <- W.Sum[gene, cell_line]
      df_per_gene[cell_line, "ULM"] <- ULM[gene, cell_line]
      df_per_gene[cell_line, "MLM"] <- MLM[gene, cell_line]
      df_per_gene[cell_line, "Viper"] <- Viper[gene, cell_line]
      df_per_gene[cell_line, "Consensus"] <- Consensus[gene, cell_line]
      df_per_gene[cell_line, "Expression"] <- Expression[gene, cell_line]
      df_per_gene[cell_line, "Gene_effect"] <- gene_effect_cancer[gene, cell_line]
    }
    
    df_per_gene$DepMap_ID <- rownames(df_per_gene)
    
    df_per_gene <- df_per_gene %>%
      pivot_longer(cols = c("ULM", "MLM", "Viper", "W.Mean", "W.Sum", "Consensus", "Expression"), 
                   names_to = "Method", values_to = "value")
    
    
    df_per_gene$Method <- factor(df_per_gene$Method, 
                                 levels = c("Expression", "Consensus", 
                                            "W.Sum", "W.Mean", "Viper", "MLM", "ULM"))

    #make scatterplots for correlations
    plot <- df_per_gene %>%
      ggplot(aes(x = Gene_effect, y = value)) +
      geom_point(size = 0.8, shape = 1) +
      geom_smooth(method = lm, se = FALSE, colour = "#454851", size = 0.8) +
      #geom_label_repel(aes(label = cell_line_name), size = 2.5, max.overlaps = 25) +
      stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
               method = "pearson", label.sep = "\n", size = 2.5, label.x.npc = "centre") +
      facet_wrap(vars(Method), nrow = 3, scales = "free") + 
      theme_bw(base_size = 14) +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
      labs(x = "Essentiality score", y = "Expression/Activity value", 
           title = paste0(gene, " - ", cancer_label, " - ", regulons)) 
           
    #save plots
    output_dir <- paste0("/home/cosmin/depmap_analysis/results/individual_genes/", regulons)
    ggsave(filename = paste0(output_dir, "/", gene, "_", cancer_label, "_", regulons, "_cor.png"), plot = plot, 
                  width = 22, height = 15, dpi = 1000, units = "cm")
           
  }
}