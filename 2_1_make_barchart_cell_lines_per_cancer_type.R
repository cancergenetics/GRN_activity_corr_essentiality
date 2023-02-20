#This script creates Fig 2A
#Input: processed Depmap data, i.e., filtered_annotated_sample_info.csv and filtered_gene_effect output from script prepare_input_data.R
#Output: Fig 2A showing how many cell lines with both expression and essentiality info are present in the Depmap 

library(tidyverse)
library(Cairo)

gene_effect <- read.csv("./depmap_input_data/wrangled_input_data/gene_effect_filtered.csv", row.names = 1)
cancer_types <- read.csv("./depmap_input_data/wrangled_input_data/filtered_annotated_sample_info.csv") %>%
  na.omit() %>%
  dplyr::filter(DepMap_ID %in% rownames(gene_effect)) %>%
  plyr::count("cancer")

ccle_cell_lines_plot <- cancer_types %>%
  mutate_at("cancer", .funs = toupper) %>%
  ggplot(aes(x = reorder(cancer, -freq) , y = freq)) +
  geom_bar(stat = "identity", fill = "#102542") +
  coord_flip() + 
  scale_y_continuous(limits = c(0, 53), breaks = scales::pretty_breaks(n = 6)) +
  geom_text(aes(label = freq, hjust = -0.25)) + 
  labs(x = "Cancer", y = "No. cell lines") + 
  theme_classic(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = "./plots/ccle_cell_line_per_cancer.pdf", plot = ccle_cell_lines_plot, 
       width = 10, height = 10, dpi = 800, units = "cm", device = cairo_pdf)
