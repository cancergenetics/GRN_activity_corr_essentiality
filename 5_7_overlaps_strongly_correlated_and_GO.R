#This script compares the genes found by three GRN methods (ARACNe, DoRothEA and GRNdb) with an R > 0.6 between consensus activity/expression and essentiality
#Additionally, it performs GO enrichment on these genes
#Input: Correlation scores produced by script 4_0_run_correlations.R for the three GRN methods
#Output: Upset plots of the overlaps between activity inferred based on ARACNe, DoRothEA, GRNdb and expresssion of genes found to be significantly correlated with essentiality
#Output: GO enrichment for the significant genes of each methods
#Output: Venn diagram of overlapping GO terms between the methods
#Output: Supplementary Fig S6A
#Output: Supplementary Fig S5C,D,E,F

library(tidyverse)
library(UpSetR)
library(WebGestaltR)
library(VennDiagram)

regulons_list <- c("aracne", "grndb", "dorothea")
cancer_list = c("aml", "blca", "brca", "coad",
                "gbm", "hnsc", "kirc", "luad",
                "paad", "stad")

#select |R| threshold for significant genes to plot in upset and to do go analysis on
threshold <- 0.4 #OR
threshold <- 0.6 #OR 
threshold <- 0.8

#get lists with genes where activity-essentiality correlations > chosen threshold 
for (cancer_label in cancer_list) {
  for (regulons in regulons_list) {
    if (regulons == "dorothea") {
      signif_genes <- read.csv(paste0("./results/results_matrices_per_combination/", 
                                                  regulons, "/", regulons, "_", cancer_label,
                                                  "/corr_scores_filtered_", regulons, "_", cancer_label, 
                                                  ".csv"), row.names = 1) %>%
        dplyr::filter(abs(Consensus.cor) > threshold) %>%
        dplyr::pull("hgnc_symbol")
      eval(parse(text = paste0("signif_", regulons, "_", cancer_label, " <- signif_genes")))
  
      } else { #regulons not dorothea
        signif_genes <- read.csv(paste0("./results/results_matrices_per_combination/",
                                        regulons, "/", regulons, "_", cancer_label, "_", cancer_label,
                                        "/corr_scores_filtered_", regulons, "_", cancer_label, "_", cancer_label,
                                        ".csv"), row.names = 1) %>%
          dplyr::filter(abs(Consensus.cor) > threshold) %>%
          dplyr::pull("hgnc_symbol")
        eval(parse(text = paste0("signif_", regulons, "_", cancer_label, " <- signif_genes")))
      }
  }
}

#get expression list
for (cancer_label in cancer_list) {
  regulons <- "aracne"
  signif_genes <- read.csv(paste0("./results/results_matrices_per_combination/",
                                  regulons, "/", regulons, "_", cancer_label, "_", cancer_label,
                                  "/corr_scores_filtered_", regulons, "_", cancer_label, "_", cancer_label,
                                  ".csv"), row.names = 1) %>%
    dplyr::filter(abs(Expression.cor) > threshold) %>%
    dplyr::pull("hgnc_symbol")
  eval(parse(text = paste0("signif_expr_", regulons, "_", cancer_label, " <- signif_genes")))
  
}

###COAD upset plots
coad_overlap <- list(ARACNe = signif_aracne_coad, DoRothEA = signif_dorothea_coad, GRNdb = signif_grndb_coad, 
                     Expression = signif_expr_aracne_coad)

upset(fromList(coad_overlap), order.by = "degree", 
      text.scale = 2, main.bar.color = "#102542",
      point.size = 3.5, line.size = 1.7, 
      sets.x.label = paste0("No. genes\nwith |R| > ", threshold)) 


###BRCA upset plots
brca_overlap <- list(ARACNe = signif_aracne_brca, DoRothEA = signif_dorothea_brca, GRNdb = signif_grndb_brca, 
                     Expression = signif_expr_aracne_brca)

upset(fromList(brca_overlap), order.by = "degree", 
      text.scale = 2, main.bar.color = "#102542",
      point.size = 3.5, line.size = 1.7, 
      sets.x.label = paste0("No. genes\nwith |R| > ", threshold)) + 
  scale_x_continuous(breaks = seq(0, 5, by = 1))


####GO analysis

#Create background lists for each of the GRN methods for the GO enrichment analysis
background_list_aracne <- c() #ARACNe
for (cancer_label in cancer_list) {
  regulons <- "aracne"
  signif_genes <- read.csv(paste0("./results/results_matrices_per_combination/",
                                  regulons, "/", regulons, "_", cancer_label, "_", cancer_label,
                                  "/corr_scores_filtered_", regulons, "_", cancer_label, "_", cancer_label,
                                  ".csv"), row.names = 1) %>% 
    dplyr::pull("hgnc_symbol")
  background_list_aracne <- c(background_list_aracne, signif_genes)
}
background_list_aracne <- unique(background_list_aracne)

background_list_dorothea <- c() #DoRothEA
for (cancer_label in cancer_list) {
  regulons <- "dorothea"
  signif_genes <- read.csv(paste0("./results/results_matrices_per_combination/", 
                                  regulons, "/", regulons, "_", cancer_label,
                                  "/corr_scores_filtered_", regulons, "_", cancer_label, 
                                  ".csv"), row.names = 1) %>%
    dplyr::pull("hgnc_symbol")
  background_list_dorothea <- c(background_list_dorothea, signif_genes)
}
background_list_dorothea <- unique(background_list_dorothea)

background_list_grndb <- c() #GRNdb
for (cancer_label in cancer_list) {
  regulons <- "grndb"
  signif_genes <- read.csv(paste0("./results/results_matrices_per_combination/",
                                  regulons, "/", regulons, "_", cancer_label, "_", cancer_label,
                                  "/corr_scores_filtered_", regulons, "_", cancer_label, "_", cancer_label,
                                  ".csv"), row.names = 1) %>% 
    dplyr::pull("hgnc_symbol")
  background_list_grndb <- c(background_list_grndb, signif_genes)
}
background_list_grndb <- unique(background_list_grndb)


#merge genes found significant in all cancer types
signif_genes_aracne <- c(signif_aracne_aml, signif_aracne_blca, signif_aracne_brca, signif_aracne_coad,
                         signif_aracne_gbm, signif_aracne_hnsc, signif_aracne_kirc, signif_aracne_luad,
                         signif_aracne_paad, signif_aracne_stad) %>%
  unique()

signif_genes_grndb <- c(signif_grndb_aml, signif_grndb_blca, signif_grndb_brca, signif_grndb_coad,
                        signif_grndb_gbm, signif_grndb_hnsc, signif_grndb_kirc, signif_grndb_luad,
                        signif_grndb_paad, signif_grndb_stad) %>%
  unique()

signif_genes_dorothea <- c(signif_dorothea_aml, signif_dorothea_blca, signif_dorothea_brca, signif_dorothea_coad,
                           signif_dorothea_gbm, signif_dorothea_hnsc, signif_dorothea_kirc, signif_dorothea_luad,
                           signif_dorothea_paad, signif_dorothea_stad) %>%
  unique()

signif_genes_expression <- c(signif_expr_aracne_aml, signif_expr_aracne_blca, signif_expr_aracne_brca, signif_expr_aracne_coad,
                             signif_expr_aracne_gbm, signif_expr_aracne_hnsc, signif_expr_aracne_kirc, signif_expr_aracne_luad,
                             signif_expr_aracne_paad, signif_expr_aracne_stad) %>%
  unique()


#GO analysis
aracne_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                                     "geneontology_Cellular_Component_noRedundant",
                                                     "geneontology_Molecular_Function_noRedundant"),
                                  interestGeneType = "genesymbol",
                                  referenceGeneType = "genesymbol",
                                  interestGene = signif_genes_aracne,
                                  referenceGene = background_list_aracne,
                                  sigMethod = "fdr", fdrThr = 0.1)

grndb_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                            "geneontology_Cellular_Component_noRedundant",
                                            "geneontology_Molecular_Function_noRedundant"),
                         interestGeneType = "genesymbol",
                         referenceGeneType = "genesymbol",
                         interestGene = signif_genes_grndb,
                         referenceGene = background_list_grndb,
                         sigMethod = "fdr", fdrThr = 0.1)

dorothea_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                           "geneontology_Cellular_Component_noRedundant",
                                           "geneontology_Molecular_Function_noRedundant"),
                        interestGeneType = "genesymbol",
                        referenceGeneType = "genesymbol",
                        interestGene = signif_genes_dorothea,
                        referenceGene = background_list_dorothea,
                        sigMethod = "fdr", fdrThr = 0.1)

expression_go <- WebGestaltR(enrichDatabase = c("geneontology_Biological_Process_noRedundant",
                                              "geneontology_Cellular_Component_noRedundant",
                                              "geneontology_Molecular_Function_noRedundant"),
                           interestGeneType = "genesymbol",
                           referenceGeneType = "genesymbol",
                           interestGene = signif_genes_expression,
                           referenceGene = background_list_aracne,
                           sigMethod = "fdr", fdrThr = 0.1)

#save GO results
write.csv(aracne_go, file = "./results/go_enrichment/aracne_consensus_activity_go_enriched_terms_R_above_0.6.csv", row.names = FALSE, quote = FALSE)
write.csv(grndb_go, file = "./results/go_enrichment/grndb_consensus_activity_go_enriched_terms_R_above_0.6.csv", row.names = FALSE, quote = FALSE)
write.csv(expression_go, file = "./results/go_enrichment/expression_go_enriched_terms_R_above_0.6.csv", row.names = FALSE, quote = FALSE)

#Venn of overlapping GO terms
venn.diagram(list(Expression = expression_go$geneSet, ARACNe = aracne_go$geneSet, GRNdb = grndb_go$geneSet), 
             lwd = 0, cex = 2, cat.cex = 3, print.mode = "raw", cat.pos = c(-25, 32, 0),
             alpha = c(0.5, 0.5, 0.5), fill = c("#6EB4D1", "#E1CE7A", "#D5573B"),
             "./plots/overlap.go.grn.methods.tiff")

#Look at the terms that are present in GRN methods, but not expression
setdiff(aracne_go$description, expression_go$description)
setdiff(grndb_go$description, expression_go$description)