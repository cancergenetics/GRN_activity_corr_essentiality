#This script processes the TCGA data so it can be used to update the GRNdb regulons

library(tidyverse)
library(impute)

#loading data downloaded from UCSC treehouse
treehouse.sample.info <- read_tsv("./tcga_data/clinical_TumorCompendium_v10_PolyA_2020-01-28.tsv")
treehouse_data <- as.data.frame(read_tsv("./tcga_data/TumorCompendium_v10_PolyA_hugo_log2tpm_58581genes_2019-07-25.tsv"))

#separate large datframe into each cancer type
aml <- treehouse.sample.info %>% #this needs extra filtering to only select TCGA samples
  dplyr::filter(treehouse.sample.info$disease == "acute myeloid leukemia" & treehouse.sample.info$site_id == "TCGA") %>%
  pull("th_sampleid")

gbm <- treehouse.sample.info %>%
  dplyr::filter(treehouse.sample.info$disease == "glioblastoma multiforme") %>%
  pull("th_sampleid")

luad <- treehouse.sample.info %>%
  dplyr::filter(treehouse.sample.info$disease == "lung adenocarcinoma") %>%
  pull("th_sampleid")

coad <- treehouse.sample.info %>%
  dplyr::filter(treehouse.sample.info$disease == "colon adenocarcinoma") %>%
  pull("th_sampleid")

brca <- treehouse.sample.info %>%
  dplyr::filter(treehouse.sample.info$disease == "breast invasive carcinoma") %>%
  pull("th_sampleid")

paad <- treehouse.sample.info %>%
  dplyr::filter(treehouse.sample.info$disease == "pancreatic adenocarcinoma") %>%
  pull("th_sampleid")

hnsc <- treehouse.sample.info %>%
  dplyr::filter(treehouse.sample.info$disease == "head & neck squamous cell carcinoma") %>%
  pull("th_sampleid")

stad <- treehouse.sample.info %>%
  dplyr::filter(treehouse.sample.info$disease == "stomach adenocarcinoma") %>%
  pull("th_sampleid")

blca <- treehouse.sample.info %>%
  dplyr::filter(treehouse.sample.info$disease == "bladder urothelial carcinoma") %>%
  pull("th_sampleid")

kirc <- treehouse.sample.info %>%
  dplyr::filter(treehouse.sample.info$disease == "kidney clear cell carcinoma") %>%
  pull("th_sampleid")


samples_list <- c("aml", "gbm", "luad", "coad", 
                  "brca", "paad", "hnsc", 
                  "stad", "blca", "kirc")

for (samples in samples_list) { 
  rna_tpm <- treehouse_data %>%
    subset(select = c("Gene", eval(parse(text = samples))))
  rownames(rna_tpm) <- rna_tpm$Gene 
  rna_tpm <- rna_tpm %>%
    subset(select = -Gene)
  
  rna_tpm[rna_tpm == 0] <- NA
  rna_tpm <- rna_tpm[which(rowMeans(!is.na(rna_tpm)) > 0.5), ] #delete rows with more than half of values = NA
  rna_tpm.mat <- as.matrix(rna_tpm)
  
  ###imputation part is inspired from the NETBiD workflow ---- from https://jyyulab.github.io/NetBID/
  #number of NAs
  sample_na_count <- apply(rna_tpm.mat, 1, function(x){length(which(is.na(x)==TRUE))})
  gene_na_count <- apply(rna_tpm.mat,2,function(x){length(which(is.na(x)==TRUE))})

  #perform imputation 
  if(sum(sample_na_count) + sum(gene_na_count)>0) rna_tpm.mat <- impute.knn(rna_tpm.mat)$data
  
  #remove lowly expressed genes
  rna_tpm.mat <- rna_tpm.mat[apply(rna_tpm.mat<=2, 1, sum)<=ncol(rna_tpm.mat)*0.90, ]
  rna_tpm_filtered <- as.data.frame(rna_tpm.mat)
  
  #save a separate matrix of the transcriptome of each cancer type 
  write.csv(rna_tpm_filtered, file = paste0("./tcga_data/", samples, "_tcga_processed.csv")) 
  print(paste("DONE", samples, "RNA-Seq"))
  }
