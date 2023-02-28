#This script produces Fig S3
#Input: for each regulons source and cancer type, the inferred activity and mRNA expression
#Output: boxplot comparing the per gene variance in activity/expression acorss all regulon sources and cancer types

library(matrixStats)
library(tidyverse)
library(Cairo)

#function to calculate per gene variance 
calculate_variance_per_gene <- function(results, method) {
  variance <- results %>%
    as.matrix() %>%
    rowVars(useNames = TRUE) %>%
    as.data.frame() %>%
    rename_at(1, ~method) %>%
    rownames_to_column("Gene")
  
  return(variance)
}

cancer_list = c("aml", "blca", "brca", "coad",
                "gbm", "hnsc", "kirc", "luad",
                "paad", "stad")

regulons_list <- c("aracne", "grndb", "dorothea")

merged_plot_variance_df <- as.data.frame(matrix(nrow = 0, ncol = 3))
colnames(merged_plot_variance_df) <- c("Gene", "method", "variance")

for (regulons_types in regulons_list) {
  for (cancer_label in cancer_list) {
    
    regulon <- paste0(regulons_types, "_", cancer_label)
    
    if (regulons_types == "dorothea") {
      regulon <- regulons_types
    } else {
      regulon <- paste0(regulons_types, "_", cancer_label)
    }
    
    input_dir <- paste0("./results/results_matrices_per_combination/",
                        regulons_types, "/", regulon, "_", cancer_label)
    
    results.W.Mean <- read.csv(paste0(input_dir, "/", regulon, "_", cancer_label, "_wmean.csv"), row.names = 1)
    results.W.Sum <- read.csv(paste0(input_dir, "/", regulon, "_", cancer_label, "_wsum.csv"), row.names = 1)
    results.ULM <- read.csv(paste0(input_dir, "/", regulon, "_", cancer_label, "_ulm.csv"), row.names = 1)
    results.MLM <- read.csv(paste0(input_dir, "/", regulon, "_", cancer_label, "_mlm.csv"), row.names = 1)
    results.Viper <- read.csv(paste0(input_dir, "/", regulon, "_", cancer_label, "_viper.csv"), row.names = 1)
    results.Consensus <- read.csv(paste0(input_dir, "/", regulon, "_", cancer_label, "_consensus.csv"), row.names = 1)
    t_gene_exp_tmp_subset <- read.csv(paste0(input_dir, "/", regulon, "_", cancer_label, "_expression.csv"), row.names = 1, check.names = FALSE)
    
    variance_mlm_per_gene <- calculate_variance_per_gene(results.MLM, "MLM")
    variance_ulm_per_gene <- calculate_variance_per_gene(results.ULM, "ULM")
    variance_w.mean_per_gene <- calculate_variance_per_gene(results.W.Mean, "W.Mean")
    variance_w.sum_per_gene <- calculate_variance_per_gene(results.W.Sum, "W.Sum")
    variance_viper_per_gene <- calculate_variance_per_gene(results.Viper, "Viper")
    variance_consensus_per_gene <- calculate_variance_per_gene(results.Consensus, "Consensus")
    variance_expression_per_gene <- calculate_variance_per_gene(t_gene_exp_tmp_subset, "Expression")
    
    variance_df <- merge(variance_expression_per_gene, variance_ulm_per_gene, by = "Gene") %>%
      merge(variance_mlm_per_gene, by = "Gene") %>%
      merge(variance_viper_per_gene, by = "Gene") %>%
      merge(variance_w.mean_per_gene, by = "Gene") %>% 
      merge(variance_w.sum_per_gene, by = "Gene") %>%
      merge(variance_consensus_per_gene, by = "Gene")
    
    plot_variance_df <- variance_df %>%
      pivot_longer(cols = c( "Expression", "MLM", "ULM", "Viper", "W.Mean", "W.Sum", "Consensus"), 
                   values_to = "variance", names_to = "method")
    
    #log2 to make plotting easier, otherwise everything will be skewed
    plot_variance_df$variance <- log2(plot_variance_df$variance)
    
    #You can uncomment these lines if you wish to produce a variance plot for each cancer type+regulon pair
    #otherwise only the aggregated plot will be produced, after the for loop
    
    #variance_plot <- ggplot(plot_variance_df, aes(x = method, y = variance)) +
     # geom_boxplot(fill = "#61a5c2", size = 0.6) +
      #theme_classic(base_size = 16) + 
      #scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
      #theme(legend.position = "none") +
      #coord_flip() +
      #labs(x = "Method", y = "log2 variance per gene", title = paste(cancer_label, regulons_types))
    
    #ggsave(filename = paste0("./plots/method_variance_per_gene_", cancer_label, "_", regulons_types, ".pdf"), 
     #      plot = variance_plot, width = 20, height = 10, dpi = 1000, units = "cm", device = cairo_pdf) 
    
    merged_plot_variance_df <- rbind(merged_plot_variance_df, plot_variance_df)
  }
}

#shows the values
plot_variance_df %>%
  group_by(method) %>%
  summarise(
    n = n(),
    mean = mean(variance),
    median = median(variance),
    sd = sd(variance))

#Viper to VIPER
merged_plot_variance_df <- merged_plot_variance_df %>% 
  mutate(method = str_replace(method, "Viper", "VIPER"))

#for ordering the method labels
merged_plot_variance_df$method <- factor(merged_plot_variance_df$method, 
                                         levels=c("Consensus", "W.Sum", "W.Mean", "VIPER", "MLM", "ULM", "Expression"))

#create plot
variance_plot_overall <- ggplot(merged_plot_variance_df, aes(x = method, y = variance)) +
  geom_boxplot(fill = "#61a5c2", size = 0.6, outlier.size = 1) +
  theme_classic(base_size = 16) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
  theme(legend.position = "none") + 
  coord_flip() +
  labs(x = "Method", y = "Log2 variance per gene", title = "Across all regulons")

#save plot
ggsave(filename = "./variance_per_gene_activity_vs_expression.pdf", plot = variance_plot_overall, 
       width = 15, height = 10, dpi = 1000, units = "cm", device = cairo_pdf) 
