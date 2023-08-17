#This script prepares gene expression data from the CCLE to construct regulons
#Input: CCLE gene expression matrix, annotated sample info data, after filtering step
#Output: gene expression for LUAD, COAD, BRCA, PAAD and their respective TF lists

library(tidyverse)

#loading raw downloaded data
gene_exp <- read.csv("./depmap_input_data/CCLE_expression.csv") %>%
  dplyr::rename(DepMap_ID = X)
filtered_annotated_sample_info <- read.csv("./depmap_input_data/wrangled_input_data/filtered_annotated_sample_info.csv")


colnames(gene_exp) <- sub("\\..*", "", colnames(gene_exp))

#remove duplicate columns
gene_exp <- gene_exp[, !duplicated(colnames(gene_exp))] 


#get a gene expression matrix for each of the cancer types used
gene_exp_luad <- gene_exp %>%
  dplyr::filter(DepMap_ID %in% filtered_annotated_sample_info[filtered_annotated_sample_info$cancer == "luad", "DepMap_ID"])
rownames(gene_exp_luad) <- gene_exp_luad$DepMap_ID
gene_exp_luad <- gene_exp_luad %>%
  dplyr::select(-DepMap_ID) %>%
  t() %>%
  as.data.frame()
gene_exp_luad$gene <- rownames(gene_exp_luad)
gene_exp_luad <- gene_exp_luad %>%
  select(gene, everything())

gene_exp_brca <- gene_exp %>%
  dplyr::filter(DepMap_ID %in% filtered_annotated_sample_info[filtered_annotated_sample_info$cancer == "brca", "DepMap_ID"])
rownames(gene_exp_brca) <- gene_exp_brca$DepMap_ID
gene_exp_brca <- gene_exp_brca %>%
  dplyr::select(-DepMap_ID) %>%
  t() %>%
  as.data.frame()
gene_exp_brca$gene <- rownames(gene_exp_brca)
gene_exp_brca <- gene_exp_brca %>%
  select(gene, everything())

gene_exp_coad <- gene_exp %>%
  dplyr::filter(DepMap_ID %in% filtered_annotated_sample_info[filtered_annotated_sample_info$cancer == "coad", "DepMap_ID"])
rownames(gene_exp_coad) <- gene_exp_coad$DepMap_ID
gene_exp_coad <- gene_exp_coad %>%
  dplyr::select(-DepMap_ID) %>%
  t() %>%
  as.data.frame()
gene_exp_coad$gene <- rownames(gene_exp_coad)
gene_exp_coad <- gene_exp_coad %>%
  select(gene, everything())

gene_exp_paad <- gene_exp %>%
  dplyr::filter(DepMap_ID %in% filtered_annotated_sample_info[filtered_annotated_sample_info$cancer == "paad", "DepMap_ID"])
rownames(gene_exp_paad) <- gene_exp_paad$DepMap_ID
gene_exp_paad <- gene_exp_paad %>%
  dplyr::select(-DepMap_ID) %>%
  t() %>%
  as.data.frame()
gene_exp_paad$gene <- rownames(gene_exp_paad)
gene_exp_paad <- gene_exp_paad %>%
  select(gene, everything())

#save results to run with ARACNe-AP
write.table(gene_exp_luad, file = "./ccle_regulons_input/luad_expr.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_exp_brca, file = "./ccle_regulons_input/brca_expr.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_exp_coad, file = "./ccle_regulons_input/coad_expr.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_exp_paad, file = "./ccle_regulons_input/paad_expr.txt", sep = "\t", quote = FALSE, row.names = FALSE)

aracne_luad <- read.csv("./regulons/aracne_regulons/luad_aracne_regulon.csv")
aracne_brca <- read.csv("./regulons/aracne_regulons/brca_aracne_regulon.csv")
aracne_coad <- read.csv("./regulons/aracne_regulons/coad_aracne_regulon.csv")
aracne_paad <- read.csv("./regulons/aracne_regulons/paad_aracne_regulon.csv")

#get a list of TFs to run ARACNe-AP
luad_tfs <- as.data.frame(unique(aracne_luad$tf))
brca_tfs <- as.data.frame(unique(aracne_brca$tf))
coad_tfs <- as.data.frame(unique(aracne_coad$tf))
paad_tfs <- as.data.frame(unique(aracne_paad$tf))

write.table(luad_tfs, file = "./ccle_regulons_input/tfs_luad.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(brca_tfs, file = "./ccle_regulons_input/tfs_brca.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(coad_tfs, file = "./ccle_regulons_input/tfs_coad.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(paad_tfs, file = "./ccle_regulons_input/tfs_paad.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)