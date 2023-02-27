#This script creates Fig. 2B,C,D
#Input: correlations between activity/expression and sensitivity
#Ouput: a csv file containing the average correlation per cancer type, for each method
#Output: The plot in Fig 2B,C,D

library(tidyverse)
library(stringr)

#select regulons
regulons <- "aracne" #OR
regulons <- "dorothea" #OR
regulons <- "grndb"

input_dir <- paste0("./results/correlation_matrices/", regulons)

cancer_list = c("aml", "blca", "brca", "coad",
                "gbm", "hnsc", "kirc", "luad",
                "paad", "stad")

methods <- c("Expression", "Consensus", "ULM", "MLM", "VIPER", "W.Mean", "W.Sum")

#load correlations matrices
VIPER <- read.csv(paste0(input_dir, "/viper_corr_pearson_", regulons, ".csv"), row.names = 1)
ULM <- read.csv(paste0(input_dir, "/ulm_corr_pearson_", regulons, ".csv"), row.names = 1)
MLM <- read.csv(paste0(input_dir, "/mlm_corr_pearson_", regulons, ".csv"), row.names = 1)
W.Mean <- read.csv(paste0(input_dir, "/wmean_corr_pearson_", regulons, ".csv"), row.names = 1)
W.Sum <- read.csv(paste0(input_dir, "/wsum_corr_pearson_", regulons, ".csv"), row.names = 1)
Consensus <- read.csv(paste0(input_dir, "/consensus_corr_pearson_", regulons, ".csv"), row.names = 1)
Expression <- read.csv(paste0(input_dir, "/expr_corr_pearson_", regulons, ".csv"), row.names = 1)

#create empty df to populate with results
cancer_method_matched_regulon <- data.frame(matrix(ncol = length(cancer_list), nrow = length(methods)))
colnames(cancer_method_matched_regulon) <- cancer_list
rownames(cancer_method_matched_regulon) <- methods

#make a column with the methods
cancer_method_matched_regulon$Method <- rownames(cancer_method_matched_regulon)

#very important: order in cancer_list must be the same as order in the activity files: VIPER, ULM, etc
for (method in methods) {
  for (i in 1:10) {
    if (regulons != "dorothea") {
      cancer_method_matched_regulon[method, i] <- eval(parse(text = method))[i,i]
    } else {
      cancer_method_matched_regulon[method, i] <- eval(parse(text = method))[i,1]
    }
  }
}

#pivot results to get them in a tidy format
cancer_method_matched_regulon_pivoted <- cancer_method_matched_regulon %>%
  pivot_longer(cols = -Method, names_to = "Cancer", values_to = "Pearsons_R")

cancer_method_matched_regulon_pivoted$Cancer <- factor(cancer_method_matched_regulon_pivoted$Cancer,
                                                        levels = cancer_list) #make factors for ggplot

cancer_method_matched_regulon_pivoted$Method <- factor(cancer_method_matched_regulon_pivoted$Method,
                                                       levels = methods) #make factors for ggplot

#separate activity methods from expression with a new column
cancer_method_matched_regulon_pivoted <- cancer_method_matched_regulon_pivoted %>%
  mutate(method_type = case_when(Method == "Expression" ~ "Expression",
                                 Method != "Expression" ~ "Activity"))

#create column that stores the number of sometimes-essential genes for which correlations were calculated
cancer_method_matched_regulon_pivoted$no_genes <- NA

#reading files depending on regulons selected
for (cancer_type in cancer_list) {
  if(regulons != "dorothea") {
    no_genes <- read.csv(paste0("./results/results_matrices_per_combination/", 
                              regulons, "/", regulons, "_", cancer_type, "_", cancer_type,
                              "/corr_scores_filtered_", regulons, "_", cancer_type,  
                              "_", cancer_type, ".csv"), row.names = 1) %>%
      nrow()
  } else {
    no_genes <- read.csv(paste0("./results/results_matrices_per_combination/", 
                                regulons, "/", regulons, "_", cancer_type,
                                "/corr_scores_filtered_", regulons, "_", cancer_type,
                                ".csv"), row.names = 1) %>%
      nrow()    
  }
  
    cancer_method_matched_regulon_pivoted[cancer_method_matched_regulon_pivoted$Cancer == cancer_type, "no_genes"] <- no_genes

    }

#saving summarised information for each regulon source 
write.csv(cancer_method_matched_regulon_pivoted, file = paste0("./results/plots_data/average_corr_per_cancer_type_", regulons, ".csv"), row.names = FALSE)

#merge strings: cancer label + number of genes used --> to appear in plot titles
cancer_method_matched_regulon_pivoted$Cancer <- 
  str_c(cancer_method_matched_regulon_pivoted$Cancer, 
        "\nN=",
        cancer_method_matched_regulon_pivoted$no_genes)

#create skittles plot
plot <- cancer_method_matched_regulon_pivoted %>% 
  mutate_at("Cancer", .funs = toupper) %>%
  ggplot(aes(x = Cancer, y = Pearsons_R, xlab)) + 
  geom_point(
    aes(colour = factor(Method),
        shape = factor(Method)), 
    size = 2,
    position = position_jitter(
      seed = 1, width = .2
    )) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8), limits = c(0.120, 0.278)) +    # set limits accordingly
  scale_colour_manual(values = c("#000000", "#009E73", "#D55E00", "#E69F00", "#56B4E9", "#CC79A7", "#0072B2"), 
                      name = "Method") +
  scale_shape_manual(values = c(17, rep(16, 6)), name = "Method") +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_line(size = 12),
        axis.ticks.x = element_blank()) +
  labs(title = paste0("Skittles plot - ", regulons, " - expr included - Pearson"), 
       y = as.expression(bquote(bar('|R|')), 
                         subtitle = paste0(regulons, " Regulons")), 
       caption = "text") 

#saving plot
ggsave(paste0("./plots/skittles_plots/", regulons, "_expr_included_skittles_plot_methods_comparison.pdf"), 
       plot = plot, width = 22.4, height = 15, 
       dpi = 1000, units = "cm", device = cairo_pdf) 