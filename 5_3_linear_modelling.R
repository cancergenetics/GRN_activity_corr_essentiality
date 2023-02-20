#This script computes linear models that predict the correlation between essentiality and activity
#It calculates the variance % (R-squared) each term contributes to the model
#Input: Average correlations between activity and gene sensitivity for each regulon type
#Input: Variance in RNA-Seq
#Output: Plot showing the variance (R-squared value) each term contributes to predicting the correlation between essentiality and activity 

library(Cairo)
library(tidyverse)
library(stringr)

#load cancer list, methods used and regulons sources
cancer_list = c("aml", "blca", "brca", "coad",
                "gbm", "hnsc", "kirc", "luad",
                "paad", "stad")

methods <- c("Expression", "Consensus", "ULM", "MLM", "VIPER", "W.Mean", "W.Sum")

regulons_types <- c("grndb")

lm <- NULL

for (regulons in regulons_types) {
  
  input_dir <- paste0("./results/correlation_matrices/", regulons)
  
  VIPER <- read.csv(paste0(input_dir, "/viper_corr_pearson_", regulons, ".csv"), row.names = 1)
  ULM <- read.csv(paste0(input_dir, "/ulm_corr_pearson_", regulons, ".csv"), row.names = 1)
  MLM <- read.csv(paste0(input_dir, "/mlm_corr_pearson_", regulons, ".csv"), row.names = 1)
  W.Mean <- read.csv(paste0(input_dir, "/wmean_corr_pearson_", regulons, ".csv"), row.names = 1)
  W.Sum <- read.csv(paste0(input_dir, "/wsum_corr_pearson_", regulons, ".csv"), row.names = 1)
  Consensus <- read.csv(paste0(input_dir, "/consensus_corr_pearson_", regulons, ".csv"), row.names = 1)
  Expression <- read.csv(paste0(input_dir, "/expr_corr_pearson_", regulons, ".csv"), row.names = 1)
  
  cancer_method_matched_regulon <- data.frame(matrix(ncol = length(cancer_list), nrow = length(methods)))
  colnames(cancer_method_matched_regulon) <- cancer_list
  rownames(cancer_method_matched_regulon) <- methods
  
  cancer_method_matched_regulon$Method <- rownames(cancer_method_matched_regulon)
  
  #very important: order in cancer_list must be the same as order in the activity files: Viper, ULM, etc
  for (method in methods) {
    for (i in 1:10) {
      if (regulons != "dorothea") {
        cancer_method_matched_regulon[method, i] <- eval(parse(text = method))[i,i]
      } else {
        cancer_method_matched_regulon[method, i] <- eval(parse(text = method))[i,1]
      }
    }
  }
  
  cancer_method_matched_regulon_pivoted <- cancer_method_matched_regulon %>%
    pivot_longer(cols = -Method, names_to = "Cancer", values_to = "Pearsons_R")
  
  cancer_method_matched_regulon_pivoted$Cancer <- factor(cancer_method_matched_regulon_pivoted$Cancer,
                                                         levels = cancer_list)
  
  cancer_method_matched_regulon_pivoted$Method <- factor(cancer_method_matched_regulon_pivoted$Method,
                                                         levels = methods)
  
  cancer_method_matched_regulon_pivoted <- cancer_method_matched_regulon_pivoted %>%
    mutate(method_type = case_when(Method == "Expression" ~ "Expression",
                                   Method != "Expression" ~ "Activity"))
  
  cancer_method_matched_regulon_pivoted$regulons_source <- regulons
  
  lm <- rbind(lm, cancer_method_matched_regulon_pivoted)
}


#compare expression vs consensus (median and p-value of Wilcoxon test)
expr_vs_consensus <- lm %>%
  dplyr::filter(Method == "Consensus" | Method == "Expression") %>%
  dplyr::filter(regulons_source == "aracne") %>%
  dplyr::select(c("Cancer", "Method", "Pearsons_R"))

#wilcoxon paired test to check if there is a significant difference between expressiona and consensus
wilcox.test(unlist(expr_vs_consensus[expr_vs_consensus$Method == "Expression", "Pearsons_R"]), 
            unlist(expr_vs_consensus[expr_vs_consensus$Method == "Consensus", "Pearsons_R"]), 
            paired = TRUE)

#look at medians to see whether expression or consensus is larger
median(unlist(expr_vs_consensus[expr_vs_consensus$Method == "Expression", "Pearsons_R"]))
median(unlist(expr_vs_consensus[expr_vs_consensus$Method == "Consensus", "Pearsons_R"]))


#update df_lm to not contain expression for the linear modelling
lm <- lm %>%
  dplyr::filter(Method != "Expression")

mean_variance <- read.csv("./results/rna_seq_mean_variance.csv")

gene_effect <- read.csv("./depmap_input_data/wrangled_input_data/gene_effect_filtered.csv", row.names = 1)
cancer_types <- read.csv("./depmap_input_data/wrangled_input_data/filtered_annotated_sample_info.csv") %>%
  na.omit() %>%
  filter(DepMap_ID %in% rownames(gene_effect)) %>%
  plyr::count("cancer")

df_lm <- merge(lm, cancer_types, by.x = "Cancer", by.y = "cancer") %>% 
  merge(mean_variance, by.x = "Cancer", by.y = "cancer_type")


#some histograms to visualize the data
hist(df_lm$Pearsons_R)
hist(df_lm$variance_per_gene)
hist(df_lm$freq)

#Compute possible models with combinations of the five terms
model1 <- lm(Pearsons_R ~ Cancer + freq + variance_per_gene + regulons_source + Method,
             data = df_lm)
model2 <- lm(Pearsons_R ~ Cancer,
             data = df_lm)
model3 <- lm(Pearsons_R ~ freq,
             data = df_lm)
model4 <- lm(Pearsons_R ~ variance_per_gene,
             data = df_lm)
model5 <- lm(Pearsons_R ~ regulons_source,
             data = df_lm)
model6 <- lm(Pearsons_R ~ Method,
             data = df_lm)
model7 <- lm(Pearsons_R ~ freq + variance_per_gene + regulons_source + Method,
             data = df_lm)
model8 <- lm(Pearsons_R ~ Cancer + freq + variance_per_gene ,
             data = df_lm)
model9 <- lm(Pearsons_R ~ freq + variance_per_gene ,
             data = df_lm)

models <- list(model1, model2, model3, model4, model5, model6, model7, model8, model9)

#renaming the models with suggestive names, containing the terms present
models_names <- c("Cancer + No. cell lines + RNA-Seq variance + Regulons type + Method", 
                  "Cancer", 
                  "No. cell lines", 
                  "RNA-Seq variance", 
                  "Regulons type", 
                  "Method", 
                  "No. cell lines + RNA-Seq variance + Regulons type + Method", 
                  "Cancer + No. cell lines + RNA-Seq variance", 
                  "No. cell lines + RNA-Seq variance")
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
ggsave("./plots/linear_model_rsq.pdf", plot = plot, 
       width = 27, height = 11, dpi = 1000, 
       units = "cm", device = cairo_pdf) 