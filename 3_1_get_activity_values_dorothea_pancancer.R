###this scripts gets VIPER, W.Mean, W.Sum, ULM, MLM and Consensus activity values at a pancancer level for each gene
###works with dorothea regulons only

library(tidyverse)
library(decoupleR)
library(binilib)
library(org.Hs.eg.db)
library(reshape2)
library(dorothea)

# loading data ----
filtered_annotated_sample_info <- read.csv("./depmap_input_data/wrangled_input_data/filtered_annotated_sample_info.csv")
gene_effect <- read.csv("./depmap_input_data/wrangled_input_data/gene_effect_filtered.csv", row.names = 1)
t_gene_exp <- read.csv("./depmap_input_data/wrangled_input_data/t_gene_exp_filtered.csv", row.names = 1, check.names = FALSE)

#only working with dorothea regulons as they are not cancer type specific
regulons_types <- "dorothea"

#functions ----
###separates results from decoupler for each method - !specify method as a string ie. "viper"
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

#load dorothea regulons ----
data("dorothea_hs", package = "dorothea")
dorothea = dorothea_hs
dorothea <- dorothea_hs %>%
  as.data.frame() %>%
  dplyr::filter(dorothea$tf %in% rownames(t_gene_exp)) %>%
  dplyr::filter(confidence %in% c("A", "B", "C")) %>%
  base::subset(select = -confidence)

regulon <- "dorothea"
cancer_label <- "pancancer"

t_gene_exp[t_gene_exp == 0] <- NA  #replace 0's with NAs
t_gene_exp <- t_gene_exp[which(rowMeans(!is.na(t_gene_exp)) > 0.8), ] #delete genes with more than half the values NAs
t_gene_exp[is.na(t_gene_exp)] = 0

t_gene_exp_mat <- as.matrix(t_gene_exp) # transform into matrix

regulon.df <- as.data.frame(eval(parse(text = regulon)))

regulon.target.df <- regulon.df %>%
  dplyr::filter(regulon.df$tf %in% row.names(t_gene_exp)) %>%
  tidyr::drop_na() #filter to keep only tfs present in expression matrix, otherwise decoupleR throws an error

print(paste("start decoupling", regulon, "+", cancer_label))

#decoupleR main chunk----

results.merged <- decouple(
  mat = t_gene_exp_mat,
  network = regulon.target.df,
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
    viper = list(.mor = "mor", pleiotropy = TRUE, verbose = FALSE)
  )
)


#separate results per each method
results.VIPER <- sep_decoupler_results(results.merged, "viper")
results.ULM <- sep_decoupler_results(results.merged, "ulm")
results.MLM <- sep_decoupler_results(results.merged, "mlm")
results.W.Sum <- sep_decoupler_results(results.merged, "wsum")
results.W.Mean <- sep_decoupler_results(results.merged, "wmean")
results.Consensus <- sep_decoupler_results(results.merged, "consensus")


common_genes <- Reduce(intersect,  #makes sure only regulatory genes in all methods are used downstream
                         list(rownames(results.VIPER),rownames(results.ULM), 
                              rownames(results.W.Sum), rownames(results.W.Mean),
                              rownames(results.MLM), rownames(results.Consensus)))

gene_effect_signif_reg_genes <- gene_effect %>% #makes sure only drivers in all methods are used downstream
  subset(select = common_genes) 


t_gene_exp_subset <- t_gene_exp %>%  #makes sure only drivers in all methods are used downstream
  dplyr::filter(rownames(t_gene_exp) %in% common_genes)


gene_effect_signif_reg_genes <- as.data.frame(t(as.matrix(gene_effect_signif_reg_genes))) #transposing

dir.create(paste0("./results/results_matrices_per_combination/", regulons_types, "/", regulon, "_", cancer_label))
output_dir <- paste0("./results/results_matrices_per_combination/", 
             regulons_types, "/", regulon, "_", cancer_label)

#save the results
write.csv(results.merged, file = paste0(output_dir, "/", regulon, "_", cancer_label, "_merged.csv"))
write.csv(results.VIPER, file = paste0(output_dir, "/", regulon, "_", cancer_label, "_viper.csv"))
write.csv(results.ULM, file = paste0(output_dir, "/", regulon, "_", cancer_label, "_ulm.csv"))
write.csv(results.MLM, file = paste0(output_dir, "/", regulon, "_", cancer_label, "_mlm.csv"))
write.csv(results.W.Sum, file = paste0(output_dir, "/", regulon, "_", cancer_label, "_wsum.csv"))
write.csv(results.W.Mean, file = paste0(output_dir, "/", regulon, "_", cancer_label, "_wmean.csv"))
write.csv(results.Consensus, file = paste0(output_dir, "/", regulon, "_", cancer_label, "_consensus.csv"))
write.csv(t_gene_exp_subset, file = paste0(output_dir, "/", regulon, "_", cancer_label, "_expression.csv"))
write.csv(gene_effect_signif_reg_genes, file = paste0(output_dir, "/", regulon, "_", cancer_label, "_gene_effect.csv"))

print(paste(regulon,"+", cancer_label, "DONE")) 
