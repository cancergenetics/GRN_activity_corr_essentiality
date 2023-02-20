#This script calculates mean correlations between activtiy and senstivity for each method/regulon source
#Input: average correlation scores per cancer type for each of the regulon sources, from the script 5_0_make_skittlesplot.R
#Output: no files are output, only values in the R console

library(tidyverse)

means_per_method <- function(input_cancer_method_matched_regulon_pivoted, method) {
  medians_and_means <- input_cancer_method_matched_regulon_pivoted %>%
    filter(method_type == "Activity") %>%
    group_by(Cancer) %>%
    summarise(median = median(Pearsons_R), mean = mean(Pearsons_R))
  
  print(medians_and_means)
  print(mean(medians_means$mean)) #global mean
  
  #mean per method
  means_per_method <- input_cancer_method_matched_regulon_pivoted %>% 
    filter(method_type == "Activity") %>%
    group_by(Method) %>%
    summarise(mean = mean(Pearsons_R))
  means_per_method$regulons <- "method"  
  return(means_per_method)
}

#reading correlations for all regulons sources
aracne_cancer_method_matched_regulon_pivoted <- read.csv("./results/plots_data/average_corr_per_cancer_type_aracne.csv")
grndb_cancer_method_matched_regulon_pivoted <- read.csv("./results/plots_data/average_corr_per_cancer_type_grndb.csv")
dorothea_cancer_method_matched_regulon_pivoted <- read.csv("./results/plots_data/average_corr_per_cancer_type_dorothea.csv")


#mean over method
means_per_method_dorothea <- dorothea_cancer_method_matched_regulon_pivoted %>% 
  filter(method_type == "Activity") %>%
  group_by(Method) %>%
  summarise(mean = mean(Pearsons_R))
means_per_method_dorothea$regulons <- "dorothea"  
mean(means_per_method_dorothea$mean)

means_per_method_aracne <- aracne_cancer_method_matched_regulon_pivoted %>% 
  filter(method_type == "Activity") %>%
  group_by(Method) %>%
  summarise(mean = mean(Pearsons_R))
means_per_method_aracne$regulons <- "aracne"
mean(means_per_method_aracne$mean)

means_per_method_grndb <- grndb_cancer_method_matched_regulon_pivoted %>% 
  filter(method_type == "Activity") %>%
  group_by(Method) %>%
  summarise(mean = mean(Pearsons_R))
means_per_method_grndb$regulons <- "grndb" 
mean(means_per_method_grndb$mean)


means_per_method <- rbind(means_per_method_aracne, means_per_method_dorothea, means_per_method_grndb) %>%
  group_by(Method) %>%
  summarise(mean_of_means = mean(mean))
