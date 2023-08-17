#This script stratifies the skittles plot by regulon size 
#Input: Correlation matrices from script 4_0_run_correlations.R
#Output: regulon sizes csvs
#Ouput: Supplementary Fig. S2

library(tidyverse)
library(dorothea)

#select regulons
regulons <- "grndb" #OR
regulons <- "dorothea"

cancer_list = c("aml", "blca", "brca", "coad",
                "gbm", "hnsc", "kirc", "luad",
                "paad", "stad")

methods <- c("Expression", "Consensus", "ULM", "MLM", "VIPER", "W.Mean", "W.Sum")

summary_by_regulon_size <- as.data.frame(matrix(nrow = 0, ncol = 9))
summary_by_regulon_size_per_gene <- as.data.frame(matrix(nrow = 0, ncol = 4))

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
      regulon <- read.csv(paste0("./regulons/", regulons, "_regulons/", cancer_label, "_", regulons, "_regulon.csv")) %>%
        dplyr::filter(tf %in% corr_driver_act_exp_matrix$hgnc_symbol)

    }
   
  #separating regulons by their size (i.e., how many target genes each regulon hub regulates)
   regulon_size <- table(regulon$tf) %>%
      as.data.frame() %>%
      mutate(Regulon_size = case_when(Freq<=20 ~ "Small regulons", 
                              Freq>20 & Freq <=100 ~ "Medium regulons", 
                              Freq>100 ~ "Large regulons")) %>%
      rename(Gene = Var1)
    
    write.csv(regulon_size, file = paste0("./regulons/regulon_sizes/", cancer_label, "_", regulons, "_regulon_size.csv"), quote = FALSE, row.names = FALSE)
    
    regulon_size$cancer <- cancer_label
    summary_by_regulon_size_per_gene <- rbind(summary_by_regulon_size_per_gene, regulon_size)
    
    #histograms of regulon targets
    plot <- ggplot(regulon_size, aes(x = Freq)) +
      geom_histogram(binwidth = 40, colour = "black", fill = "#93C48B") + 
      theme_bw(base_size = 20) + 
      labs(x = "Regulon size", Title = cancer_label)
    
    #save histograms
    ggsave(paste0("./plots/regulon_size_", cancer_label, "_", regulons, ".png"), plot = plot, 
           width = 20, height = 12, dpi = 1000, 
           units = "cm", ) 
    
    corr_driver_act_exp_matrix <- corr_driver_act_exp_matrix %>%
      dplyr::left_join(regulon_size, by = c("hgnc_symbol" = "Gene")) %>%
      dplyr::mutate(across(where(is.numeric), abs)) %>%
      dplyr::select(c("ULM.cor", "MLM.cor", "VIPER.cor", "W.Mean.cor", "W.Sum.cor", "Consensus.cor", "Expression.cor", "Regulon_size")) %>%
      dplyr::summarise(across(everything(), mean), .by = Regulon_size)
    
    corr_driver_act_exp_matrix$cancer_type <- cancer_label
    summary_by_regulon_size <- rbind(summary_by_regulon_size, corr_driver_act_exp_matrix)
    
    pivoted_summary_by_regulon_size <- summary_by_regulon_size %>%
      rename(c("ULM" = "ULM.cor", "MLM" = "MLM.cor", "VIPER" = "VIPER.cor", "W.Mean" = "W.Mean.cor",
               "W.Sum" = "W.Sum.cor", "Consensus" = "Consensus.cor", "Expression" = "Expression.cor")) %>%
      pivot_longer(cols = c("ULM", "MLM", "VIPER", "W.Mean", "W.Sum", "Consensus", "Expression"),
                   values_to = "Pearsons_R", names_to = "Method")
    
    pivoted_summary_by_regulon_size$Cancer <- factor(pivoted_summary_by_regulon_size$cancer_type,
                                                     levels = cancer_list) #make factors for ggplot
    
    pivoted_summary_by_regulon_size$Method <- factor(pivoted_summary_by_regulon_size$Method,
                                                     levels = methods) #make factors for ggplot
    
    #producing stratified skittles plots
    plot <- pivoted_summary_by_regulon_size %>%
      mutate_at("cancer_type", .funs = toupper) %>%
      ggplot(aes(x = cancer_type, y = Pearsons_R, xlab)) +
      geom_point(aes(colour = factor(Method),
                     shape = factor(Method)),
                 size = 2,
                 position = position_jitter(
                   seed = 1, width = .2)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 8), limits = c(0.085, 0.4)) +    # set limits accordingly
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
      facet_wrap(~Regulon_size, nrow = 3)
    
    ggsave(paste0("./plots/skittles_plots/", regulons, "_by_regulon_size.pdf"),
           plot = plot, width = 22.4, height = 28,
           dpi = 1000, units = "cm", device = cairo_pdf)
    }

#save outputs for linear modelling
write.csv(pivoted_summary_by_regulon_size, file = paste0("./regulons/regulon_sizes/summary_", regulons, "_regulons_size.csv"), quote = FALSE, row.names = FALSE)
write.csv(summary_by_regulon_size_per_gene, file = paste0("./regulons/regulon_sizes/summary_per_gene", regulons, "_regulons_size.csv"), quote = FALSE, row.names = FALSE)
