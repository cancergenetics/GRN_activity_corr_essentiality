#this script modifies the data downloaded from DepMap for downstream use
#it removes entrez ids, duplicate entries and only keeps cell lines present in both expression and essentiality datasets
#it creates a label for every cancer type we analyze

#Input: 21Q4 Depmap uprocessed data downloaded at the start of the pipeline, i.e., gene expression, gene effect and cell line information
#Output: processed Depmap data for downstream use

library(tidyverse)

#loading raw downloaded data
gene_exp <- read.csv("./depmap_input_data/CCLE_expression.csv") %>%
  dplyr::rename(DepMap_ID = X)
sample.info <- read.csv("./depmap_input_data/sample_info.csv") 
gene_effect <- read.csv("./depmap_input_data/CRISPR_gene_effect.csv")


#Data wrangling

#remove entrez ids
colnames(gene_effect) <- sub("\\..*", "", colnames(gene_effect))
colnames(gene_exp) <- sub("\\..*", "", colnames(gene_exp))

#remove duplicate columns
gene_exp <- gene_exp[, !duplicated(colnames(gene_exp))] 
gene_effect <- gene_effect[, !duplicated(colnames(gene_effect))] 

#keep cell lines present in both the expression and the essentiality datasets
gene_exp <- gene_exp %>%    
  dplyr::filter(gene_exp$DepMap_ID %in% gene_effect$DepMap_ID)
gene_effect <- gene_effect %>% 
  dplyr::filter(gene_effect$DepMap_ID %in% gene_exp$DepMap_ID)

t_gene_exp <- as.data.frame(t(gene_exp)) #transpose to have genes as rows, samples as columns
colnames(t_gene_exp) <- t_gene_exp["DepMap_ID", ]  #set column names as cell lines depmap ids
t_gene_exp <- t_gene_exp[-1, ] #remove column with depmap ids
t_gene_exp[] <- as.data.frame(sapply(t_gene_exp, as.numeric)) #make sure dataframe is numeric

#keep genes present in both the expression and the essentiality datasets
t_gene_exp <- t_gene_exp %>%   
  dplyr::filter(rownames(t_gene_exp) %in% colnames(gene_effect))
gene_effect <- gene_effect %>%
  dplyr::select(c("DepMap_ID", rownames(t_gene_exp)))

#get labels of cancer cell lines used
filtered_annotated_sample_info <- sample.info %>%
  mutate(cancer = case_when(
    lineage_subtype == "AML" ~ "aml",
    Subtype == "Glioblastoma" ~ "gbm",
    Subtype == "Non-Small Cell Lung Cancer (NSCLC), Adenocarcinoma" ~ "luad",
    primary_disease == "Colon/Colorectal Cancer" & Subtype == "Adenocarcinoma" ~ "coad",
    primary_disease == "Breast Cancer" ~ "brca", 
    lineage_sub_subtype == "exocrine_adenocarcinoma" ~ "paad",
    primary_disease == "Head and Neck Cancer" & startsWith(Subtype, "Squamous Cell Carcinoma") ~ "hnsc",
    lineage_subtype == "gastric_adenocarcinoma" ~ "stad",
    lineage_subtype == "bladder_carcinoma" ~ "blca", 
    lineage_subtype == "renal_cell_carcinoma" ~ "kirc"
  )) %>%
  subset(select = c("DepMap_ID", "cancer"))

#save modified files for downstream use
write.csv(filtered_annotated_sample_info, file = "./depmap_input_data/wrangled_input_data/filtered_annotated_sample_info.csv", row.names = FALSE)
write.csv(gene_effect, file = "./depmap_input_data/wrangled_input_data/gene_effect_filtered.csv", row.names = FALSE)
write.csv(t_gene_exp, file = "./depmap_input_data/wrangled_input_data/t_gene_exp_filtered.csv", row.names = TRUE)
