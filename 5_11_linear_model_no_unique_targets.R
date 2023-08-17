#This script computes linear models that predict the correlation between essentiality and activity
#It calculates the variance % (R-squared) each term contributes to the model
#Input: Average correlations between activity and gene sensitivity for each regulon type
#Input: Variance in RNA-Seq
#Output: Plot showing the variance (R-squared value) each term contributes to predicting the correlation between essentiality and activity 
#Ouput: Supplementary Fig. S4B
library(Cairo)
library(tidyverse)
library(stringr)

#load cancer list, methods used and regulons sources
cancer_list = c("aml", "blca", "brca", "coad",
                "gbm", "hnsc", "kirc", "luad",
                "paad", "stad")

methods <- c("Expression", "Consensus", "ULM", "MLM", "VIPER", "W.Mean", "W.Sum")

regulons_types <- c("grndb", "dorothea")


lm <- NULL

for (regulons in regulons_types) {
  
  summary_unique_targets <- read.csv(paste0("./regulons/unique_targets/summary_", regulons, "_unique_targets.csv"))
  
  summary_unique_targets$Cancer <- factor(summary_unique_targets$Cancer,
                                          levels = cancer_list)
  
  summary_unique_targets$Method <- factor(summary_unique_targets$Method,
                                          levels = methods)
  
  summary_unique_targets <- summary_unique_targets %>%
    mutate(method_type = case_when(Method == "Expression" ~ "Expression",
                                   Method != "Expression" ~ "Activity"))
  
  summary_unique_targets$regulons_source <- regulons
  
  lm <- rbind(lm, summary_unique_targets)
}


#update df_lm to not contain expression for the linear modelling
lm <- lm %>%
  dplyr::filter(Method != "Expression")

cancer_types <- read.csv("./depmap_input_data/wrangled_input_data/filtered_annotated_sample_info.csv") %>%
  na.omit() %>%
  filter(DepMap_ID %in% rownames(gene_effect)) %>%
  plyr::count("cancer")

df_lm <- merge(lm, cancer_types, by.x = "Cancer", by.y = "cancer")

#some histograms to visualize the data
hist(df_lm$Pearsons_R)
hist(df_lm$variance_per_gene)
hist(df_lm$freq)

#Compute possible models with combinations of the five terms
model1 <- lm(Pearsons_R ~ Cancer + ratio + freq + regulons_source + Method,
             data = df_lm)
model2 <- lm(Pearsons_R ~ Cancer,
             data = df_lm)
model3 <- lm(Pearsons_R ~ freq,
             data = df_lm)
model4 <- lm(Pearsons_R ~ ratio,
             data = df_lm)
model5 <- lm(Pearsons_R ~ regulons_source,
             data = df_lm)
model6 <- lm(Pearsons_R ~ Method,
             data = df_lm)
model7 <- lm(Pearsons_R ~ freq + ratio + regulons_source + Method,
             data = df_lm)
model8 <- lm(Pearsons_R ~ Cancer + freq + ratio ,
             data = df_lm)
model9 <- lm(Pearsons_R ~ freq + ratio ,
             data = df_lm)

models <- list(model1, model2, model3, model4, model5, model6, model7, model8, model9)

#renaming the models with suggestive names, containing the terms present
models_names <- c("Cancer + No. cell lines + No. unique targets + Regulons type + Method", 
                  "Cancer", 
                  "No. cell lines", 
                  "No. unique targets", 
                  "Regulons type", 
                  "Method", 
                  "No. cell lines + No. unique targets + Regulons type + Method", 
                  "Cancer + No. cell lines + No. unique targets", 
                  "No. cell lines + No. unique targets")
names(models) <- models_names

#prepare dataframe for plotting the R-squared values
rsquare_plot <- as.data.frame(matrix(nrow = 9, ncol = 1))
colnames(rsquare_plot) <- c("Adjusted_R_sq")
rownames(rsquare_plot) <- models_names

#add adjusted R-sq for each model
for (model_label in models_names) {
  rsquare_plot[model_label, "Adjusted_R_sq"] <- summary(models[[model_label]])[["adj.r.squared"]]
}

rsquare_plot$model_terms <- rownames(rsquare_plot)

#multiply R-squared by 100, to get as percentage
rsquare_plot$Adjusted_R_sq <- rsquare_plot$Adjusted_R_sq * 100

#plot R-squared for each model
plot <- rsquare_plot %>%
  ggplot() +
  geom_point(aes(x = Adjusted_R_sq, y = reorder(model_terms, Adjusted_R_sq)), 
             size = 3, colour = "#1B263B") + 
  theme_bw(base_size = 14) + 
  scale_x_continuous(breaks = seq(0, 80, 10), limits = c(-1, 80)) +
  labs(x = "Adjusted R-squared (%)", y = "Model terms") +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank())

#save plot
ggsave("./plots/linear_model_rsq_unique_targets.pdf", plot = plot, 
       width = 27, height = 11, dpi = 1000, 
       units = "cm", device = cairo_pdf) 
