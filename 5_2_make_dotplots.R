library(tidyverse)
library(cowplot)
library(reshape2)
library(Cairo)

#select regulons
#only works with aracne and grndb regulons, because dorothea regulons are not cancer-type specific
regulons <- "aracne" #OR
#regulons <- "grndb"

#setting working directory to load necessary files
input_dir <- paste0("./results/correlation_matrices/", regulons)
  
#this function creates plots like the one in Fig 3  
plot_dots <- function(matrix, plotname, label, regulon_method){ #plot for activity/expression
  melted_matrix <- reshape2::melt(as.matrix(matrix))
  
  melted_matrix <- melted_matrix %>% 
    mutate(regulon_type = c('Cancer type-mismatched', 'Cancer type-matched')[1+str_detect(Var2, as.character(Var1))])
  
  melted_matrix <- melted_matrix %>%
    group_by(Var1) %>%
    mutate(rank = rank(-value))
  
  cancer_specific_ranks_regulons <- melted_matrix %>%
    filter(regulon_type == "Cancer type-matched") %>%
    pull(rank)
  
  cancer_non_specific_ranks_regulons <- melted_matrix %>%
    filter(regulon_type != "Cancer type-matched") %>%
    pull(rank)
  
  rank_test_result <- wilcox.test(cancer_specific_ranks_regulons, cancer_non_specific_ranks_regulons, 
                                  paired = FALSE, alternative =  "two.sided")
  
  plot <- melted_matrix %>%
    mutate_at("Var1", .funs = toupper) %>%
    ggplot(aes(x = Var1, y = value, xlab)) + 
    geom_point(
      aes(colour = factor(regulon_type)), 
      size = 2,
      #alpha = .2,
      position = position_jitter(
        seed = 1, width = .2
      )) +
    annotate("text", label = paste0("Unpaired Two-Samples Wilcoxon Test \n p-value = ", 
                                    round(rank_test_result[["p.value"]], digits = 3)), 
             x = 8.3, y = 0.235) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0.1, 0.24)) +
    scale_colour_manual(values = c("#F46036", "#04080F"), 
                        name = "Regulons Type", 
                        labels = c('Cancer type-\nmatched', 'Cancer type-\nmismatched')) +
    theme_bw(base_size = 16) +
    theme(panel.grid.minor.y = element_blank(), 
          panel.grid.major.x = element_line(size = 12),
          axis.ticks.x = element_blank()) +
    labs(title = paste0("Dotplot - ", label, " - ", regulon_method, " - Pearson"), 
         x = "Cancer", 
         y = as.expression(bquote(bar('|R'[.(label)]*'|')),  
                           subtitle = paste0(regulon_method, " Regulons")),
         caption = paste0("Unpaired Two-Samples Wilcoxon Test p-value = ", 
                          round(rank_test_result[["p.value"]], digits = 3))) 
  ggsave(plotname, plot = plot, width = 24, height = 15, dpi = 1000, units = "cm", device = cairo_pdf) 
  
  return(plot)
}  
  
#read files
VIPER <- read.csv(paste0(input_dir, "/viper_corr_pearson_", regulons, ".csv"), row.names = 1)
ULM <- read.csv(paste0(input_dir, "/ulm_corr_pearson_", regulons, ".csv"), row.names = 1)
MLM <- read.csv(paste0(input_dir, "/mlm_corr_pearson_", regulons, ".csv"), row.names = 1)
W.Mean <- read.csv(paste0(input_dir, "/wmean_corr_pearson_", regulons, ".csv"), row.names = 1)
W.Sum <- read.csv(paste0(input_dir, "/wsum_corr_pearson_", regulons, ".csv"), row.names = 1)
Consensus <- read.csv(paste0(input_dir, "/consensus_corr_pearson_", regulons, ".csv"), row.names = 1)
Expression <- read.csv(paste0(input_dir, "/expr_corr_pearson_", regulons, ".csv"), row.names = 1)

#vector of methods for which plots will be created
methods <- c("ULM", "MLM", "VIPER", "W.Mean", "W.Sum", "Consensus", "Expression")

for (method in methods) {
  plot_method <- plot_dots(matrix = eval(parse(text = method)),
                           plotname = paste0("./plots/dotplots/dotplot_pearson_", method, "_", regulons, ".pdf"),
                           label = method,
                           regulon_method = regulons)
}   