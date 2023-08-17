#This script stratifies the skittles plot by number of unique targets 
#Input: Correlation matrices from script 4_0_run_correlations.R
#Output: unique targets csvs
#Ouput: Supplementary Fig. S3

library(tidyverse)
library(dorothea)

#select regulons
regulons <- "grndb" #OR
regulons <- "dorothea"

cancer_list = c("aml", "blca", "brca", "coad",
                "gbm", "hnsc", "kirc", "luad",
                "paad", "stad")

methods <- c("Expression", "Consensus", "ULM", "MLM", "VIPER", "W.Mean", "W.Sum")

summary_by_unique_targets <- as.data.frame(matrix(nrow = 0, ncol = 9))
summary_unique_targets_per_gene <- as.data.frame(matrix(nrow = 0, ncol = 6))

#read files according to type of analysis chosen: pancancer or per cancer type

for (cancer_label in cancer_list) {
    if (regulons == "dorothea") {
      corr_driver_act_exp_matrix <- read.csv(paste0("./results/results_matrices_per_combination/", 
                                                    regulons, "/", regulons, "_", cancer_label,
                                                    "/corr_scores_filtered_", regulons, "_", cancer_label, 
                                                    ".csv"), row.names = 1)
      regulon <- dorothea_hs %>%
        as.data.frame() %>%
        filter(confidence %in% c("A", "B", "C")) %>%
        subset(select = -confidence) %>%
        dplyr::filter(tf %in% corr_driver_act_exp_matrix$hgnc_symbol)
      
    } else { #regulons not dorothea
      corr_driver_act_exp_matrix <- read.csv(paste0("./results/results_matrices_per_combination/", 
                                                    regulons, "/", regulons, "_", cancer_label, "_", cancer_label, 
                                                    "/corr_scores_filtered_", regulons, "_", cancer_label, "_", cancer_label,
                                                    ".csv"), row.names = 1)
      regulon <- read.csv(paste0("./regulons/", regulons, "_regulons/", cancer_label, "_", regulons, "_regulon.csv"))
      }
    
    targets_list <- unique(regulon$target)
    regulon$unique_target <- NA
    for (target_label in targets_list) {
      occurences <- regulon %>%
        dplyr::filter(target == target_label) %>%
        nrow()
      if (occurences == 1) {
        regulon[regulon[ ,"target"] == target_label, "unique_target"] <- "yes"
      } else {
        regulon[regulon[ ,"target"] == target_label, "unique_target"] <- "no"
      }
    }
    
    regulon_unique_targets_ratio <- regulon %>%
      dplyr::filter(tf %in% corr_driver_act_exp_matrix$hgnc_symbol) %>%
      dplyr::group_by(tf) %>%
      dplyr::summarise(unique_targets_count = sum(unique_target == "yes"),
                       total_targets = n()) %>%
      dplyr::mutate(ratio = unique_targets_count/total_targets * 100) %>%
      dplyr::mutate(Unique_targets = case_when(ratio == 0 ~ "No unique targets", 
                                               ratio <= 10 ~ "\u226410% of targets unique", 
                                               ratio > 10 ~ ">10% of targets unique"))
    
    write.csv(regulon_unique_targets_ratio, file = paste0("./regulons/unique_targets/", cancer_label, "_", regulons, "_unique_targets.csv"), quote = FALSE, row.names = FALSE)
    
    regulon_unique_targets_ratio$Cancer <- cancer_label
    
    corr_driver_act_exp_matrix_per_gene <- corr_driver_act_exp_matrix %>%
      dplyr::left_join(regulon_unique_targets_ratio, by = c("hgnc_symbol" = "tf")) %>%
      dplyr::mutate(across(where(is.numeric), abs)) %>%
      dplyr::select(c("hgnc_symbol", "ULM.cor", "MLM.cor", "VIPER.cor", "W.Mean.cor", 
                      "W.Sum.cor", "Consensus.cor", "Expression.cor", "Unique_targets", "ratio"))
    
    corr_driver_act_exp_matrix <- corr_driver_act_exp_matrix_per_gene %>% 
      dplyr::select(-hgnc_symbol) %>%
      dplyr::summarise(across(everything(), mean), .by = Unique_targets)
    
    corr_driver_act_exp_matrix$cancer_type <- cancer_label
    corr_driver_act_exp_matrix_per_gene$cancer_type <- cancer_label
    
    summary_by_unique_targets <- rbind(summary_by_unique_targets, corr_driver_act_exp_matrix)

    summary_unique_targets_per_gene <- rbind(summary_unique_targets_per_gene, corr_driver_act_exp_matrix_per_gene)
    
    pivoted_summary_by_unique_targets <- summary_by_unique_targets %>%
      dplyr::rename(c("ULM" = "ULM.cor", "MLM" = "MLM.cor", "VIPER" = "VIPER.cor", "W.Mean" = "W.Mean.cor",
               "W.Sum" = "W.Sum.cor", "Consensus" = "Consensus.cor", "Expression" = "Expression.cor")) %>%
      pivot_longer(cols = c("ULM", "MLM", "VIPER", "W.Mean", "W.Sum", "Consensus", "Expression"),
                   values_to = "Pearsons_R", names_to = "Method")
    
    pivoted_summary_by_unique_targets$Cancer <- factor(pivoted_summary_by_unique_targets$cancer_type,
                                                       levels = cancer_list) #make factors for ggplot
    
    pivoted_summary_by_unique_targets$Method <- factor(pivoted_summary_by_unique_targets$Method,
                                                       levels = methods) #make factors for ggplot
    
    pivoted_summary_by_unique_targets$Unique_targets <- factor(pivoted_summary_by_unique_targets$Unique_targets, 
                                                               levels = c("No unique targets", "\u226410% of targets unique", ">10% of targets unique"))
    
    #save pivoted_summary_by_unique_targets for liniar modelling
    
    plot <- pivoted_summary_by_unique_targets %>%
      mutate_at("cancer_type", .funs = toupper) %>%
      ggplot(aes(x = cancer_type, y = Pearsons_R, xlab)) +
      geom_point(aes(colour = factor(Method),
                     shape = factor(Method)),
                 size = 2,
                 position = position_jitter(
                   seed = 1, width = .2)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 8), limits = c(0, 0.76)) +    # set limits accordingly
      scale_colour_manual(values = c("#000000", "#009E73", "#D55E00", "#E69F00", "#56B4E9", "#CC79A7", "#0072B2"),
                          name = "Method") +
      scale_shape_manual(values = c(17, rep(16, 6)), name = "Method") +
      theme_bw(base_size = 16) +
      theme(panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_line(size = 12),
            axis.ticks.x = element_blank()) +
      labs(title = paste0("Skittles plot - ", regulons, " - expr included - Pearson"),
           x = "Cancer",
           y = as.expression(bquote(bar('|R|')),
                             subtitle = paste0(regulons, " Regulons")),
           caption = "text") +
      facet_wrap(~Unique_targets, nrow = 3)
    
    ggsave(paste0("./plots/skittles_plots/", regulons, "_by_unique_targets.pdf"),
           plot = plot, width = 22.4, height = 28,
           dpi = 1000, units = "cm", device = cairo_pdf)
    
    }


pivoted_summary_unique_targets_per_gene <- summary_unique_targets_per_gene %>%
  dplyr::rename(c("ULM" = "ULM.cor", "MLM" = "MLM.cor", "VIPER" = "VIPER.cor", "W.Mean" = "W.Mean.cor",
                  "W.Sum" = "W.Sum.cor", "Consensus" = "Consensus.cor", "Expression" = "Expression.cor", 
                  "Cancer" = "cancer_type")) %>%
  pivot_longer(cols = c("ULM", "MLM", "VIPER", "W.Mean", "W.Sum", "Consensus", "Expression"),
               values_to = "Pearsons_R", names_to = "Method")
  
#save outputs for linear modelling
write.csv(pivoted_summary_by_unique_targets, file = paste0("./regulons/unique_targets/summary_", regulons, "_unique_targets.csv"), quote = FALSE, row.names = FALSE)
write.csv(pivoted_summary_unique_targets_per_gene, file = paste0("./regulons/unique_targets/summary_per_gene_", regulons, "_unique_targets.csv"), quote = FALSE, row.names = FALSE)