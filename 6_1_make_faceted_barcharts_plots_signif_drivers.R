#This script creates Fig 4B,C,D,S2
#Input: .csv files output by script 6_0
#Output: barcharts counting absolute correlations 
        #barcharts showing the two types of correlations: increase in sensitivity with increase in activity/expression
                                                         #increase in sensitivity with deacrease in activity/expression

library(stringr)
library(magrittr)
library(tidyverse)
library(Cairo)

#select regulons
regulons <- "aracne" #OR
regulons <- "grndb" #OR
regulons <- "dorothea"

#select cancer list analysis
cancer_list = c("aml", "blca", "brca", "coad",
                "gbm", "hnsc", "kirc", "luad",
                "paad", "stad") #OR
cancer_list = c("pancancer") #only works with dorothea regulons

#1. creating plots for absolute correlations
barchart_absolute <- NULL

#merge results for each cancer type
for (cancer in cancer_list) {
  barchart_absolute_cancer <- read.csv(paste0("./results/barcharts_data/", 
                                     cancer, "_", regulons, "_number_signif_corr.csv")) %>%
    mutate_at("Cancer", .funs = toupper)
  barchart_absolute <- rbind(barchart_absolute, barchart_absolute_cancer)  
}

#add number of genes to each cancer as a string so it can be in each plot's header
barchart_absolute$Cancer <- str_c(barchart_absolute$Cancer,
                         '   |   No. genes = ', barchart_absolute$total_drivers)

#sorting by best performing method
vector_sorted_best_method <- barchart_absolute %>%
  group_by(Method) %>%
  summarise(sum(Percentage_signif)) %>%
  set_colnames(c("Method", "sum_percentage_signif")) %>%
  arrange(sum_percentage_signif, decreasing = TRUE) %>%
  pull(Method) %>%
  as.vector()

#make factors for ggplot to arrange in order of best to worst performing method
barchart_absolute$Method <- factor(bar_signif_correlations$Method, 
                          levels = vector_sorted_best_method)

#sets parameters for plot according to type of analysis chosen: either pancancer or per cancer type
if (cancer_list == "pancancer") {
  height <- 7.5
  width <- 15
  filename_absolute <- paste0("./plots/barcharts/1percent_", regulons, "_barcharts_5methods.pdf")
  filename_pos_vs_neg <- paste0("./plots/barcharts/1percent_", regulons, "_pos_vs_neg_5methods.pdf")
  axis_label_absolute <- "Genes significantly correlated\nwith CRISPR gene sensitivity (%)"
  axis_label_pos_vs_neg <- "Genes significantly correlated\nwith CRISPR gene sensitivity (count)"
} else {
  height <- 10
  width <- 38
  filename_absolute <- paste0("./plots/barcharts/", regulons, "_barcharts_5methods.pdf")
  filename_pos_vs_neg <- paste0("./plots/barcharts/", regulons, "_pos_vs_neg_5methods.pdf")
  axis_label_absolute <- "Genes significantly correlated with CRISPR gene sensitivity (%)"
  axis_label_pos_vs_neg <- "Genes significantly correlated with CRISPR gene sensitivity (count)"
}

#faceted barchart absolute correlations
plot_barchart_absolute <- barchart_absolute %>%
  ggplot() +
  geom_bar(aes(fill = Pearsons_R, x = Method, y = Percentage_signif),
           stat = "identity", position = position_stack(reverse = FALSE)) +
  coord_flip() +
  scale_fill_manual(values = c("#A9D6E5", "#61A5C2", "#2C7DA0", "#01497C"), 
                    name = "Pearson's R") +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  geom_text(x = 25, y = 1, label = barchart_absolute$total_drivers) +
  labs(title = paste0(regulons, " regulons"), x = "Method", y = axis_label_absolute) +
  facet_wrap(vars(Cancer), nrow = 2)  


#2. creating plots for pos vs neg correlations
barchart_pos_vs_neg <- NULL

for (cancer in cancer_list) {
  barchart_pos_vs_neg_cancer <- read.csv(paste0("./results/barcharts_data/", 
                                     cancer, "_", regulons, "_pos_vs_neg_corr.csv"))
  barchart_pos_vs_neg <- rbind(barchart_pos_vs_neg, barchart_pos_vs_neg_cancer)  
}

barchart_pos_vs_neg$count <- ifelse(barchart_pos_vs_neg$Pearsons_R == "R≤0", -1 * barchart_pos_vs_neg$count, barchart_pos_vs_neg$count)

barchart_pos_vs_neg$Method <- factor(barchart_pos_vs_neg$Method, 
                          levels = vector_sorted_best_method)

#faceted barchart pos vs neg correlations
plot_barchart_pos_vs_neg <-
  barchart_pos_vs_neg %>%
  ggplot() +
  geom_bar(aes(x = Method, y = count, fill = Pearsons_R), stat = "identity", 
           position = position_stack(reverse = FALSE)) +
  scale_fill_manual(values = c("#F7B05B", "#C2E1C2"), 
                    name = "Pearson's R", 
                    labels=c(expression("\u2191 expression\n\u2191 sensitivity"), 
                             expression("\u2191 expression\n\u2193 sensitivity"))) +
  geom_text(data = subset(barchart_pos_vs_neg, count != 0), 
            aes(x = Method, y = count, label = abs(count)), 
            hjust = ifelse(subset(barchart_pos_vs_neg, count != 0)$count > 0, 1, -0.1)) +
  scale_y_continuous(labels = abs) +
  coord_flip() +
  theme_bw(base_size = 16) + 
  theme(panel.grid = element_blank(),
        legend.text.align = 0) +
  labs(title = paste0(regulons, " regulons"), x = "Method", y = axis_label_pos_vs_neg) +
  facet_wrap(vars(Cancer), nrow = 2)  

#saving plots
ggsave(filename = filename_absolute, plot = plot_barchart_absolute, 
       width = width, height = height, 
       dpi = 1000, units = "cm", device = cairo_pdf) 

ggsave(filename = filename_pos_vs_neg, plot = plot_barchart_pos_vs_neg, 
       width = width, height = height, 
       dpi = 1000, units = "cm", device = cairo_pdf)
