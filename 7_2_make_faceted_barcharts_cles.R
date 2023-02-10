library(tidyverse)
library(magrittr)
library(Cairo)

#choose one: 
regulons <- "aracne" #OR
regulons <- "grndb" #OR
regulons <- "dorothea"

barchart <- NULL

for (cancer in cancer_list) {
  barchart_cancer <- read.csv(paste0("./results/cles/", regulons, "/",
                                     cancer, "_", regulons, "_number_signif_cles.csv"))  %>%
    mutate_at("Cancer", .funs = toupper)
  barchart <- rbind(barchart, barchart_cancer)  
  
}


vector_sorted_best_method <- barchart %>%
  group_by(Method) %>%
  summarise(sum(Percentage_signif)) %>%
  set_colnames(c("Method", "sum_percentage_signif")) %>%
  arrange(sum_percentage_signif, decreasing = TRUE) %>%
  pull(Method) %>%
  as.vector()


barchart$Method <- factor(barchart$Method, 
                          levels = vector_sorted_best_method)

#merge strings in columns so number of genes used appears in barchart title
barchart$Cancer <- str_c(barchart$Cancer, '   |   No. genes = ', barchart$total_drivers)

#sets parameters for plot according to type of analysis chosen: either pancancer or per cancer type
if (cancer_list == "pancancer") {
  height <- 7.5
  width <- 15
  filename <- paste0("./results/plots/barcharts/1percent_", regulons , "_cles_barcharts_5methods.pdf")
  axis_label <- "Conditionally essential genes\npredicted by activity/expression (%)"
} else {
  height <- 10
  width <- 38
  filename <- paste0("./results/plots/barcharts/", regulons , "_cles_barcharts_5methods.pdf")
  axis_label <- "Conditionally essential genes predicted by activity/expression (%)"
}


#create barchart
plot_barchart <- barchart %>%
  ggplot() +
  geom_bar(aes(fill = CLES, x = Method, y = Percentage_signif),
           stat = "identity", position = position_stack(reverse = FALSE)) +
  coord_flip() +
  scale_fill_manual(values = c("#ea8c55", "#ad2e24", "#540804"), 
                    name = "CLES") +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  geom_text(x = 25, y = 1, label = barchart$total_drivers) +
  labs(title = paste0(regulons, " regulons"), x = "Method", y = axis_label) +
  facet_wrap(vars(Cancer), nrow = 2)  

#save barcharts as pdf
ggsave(filename = filename, plot = plot_barchart, 
       width = width, height = height, 
       dpi = 1000, units = "cm", device = cairo_pdf)
