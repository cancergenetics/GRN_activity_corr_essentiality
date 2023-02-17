# GRN_activity_corr_essentiality

The scripts available perform the analysis presented in the manuscript by Tudose, Bond and Ryan: Gene essentiality in cancer is better predicted by mRNA abundance than by gene regulatory network-inferred activity

Each script can be run in Rstudio. However, some of the scripts require more computational time, so it is recommended the analysis is at least partly run on a server using nohup and R CMD BATCH. This is especially the case for steps 1.0-4.0

Whilst most of the code is fully automated, the user should select at the start of each run the regulons for analysis (ARACNe, DoRothEA or GRNdb) or comment out the ones not necessary.

For reproducibility, the package versions used in the original analysis are added via packrat. 


| Script                                              | Figures           | Brief description                                        |
|:----------------------------------------------------|:------------------|:---------------------------------------------------------|
| 2_1_make_barchart_cell_lines_per_cancer_type.R      | Fig. 2A           | Plot number of cell lines used for each cancer type in the CCLE |
| 5_0_make_skittlesplot.R                             | Fig. 2B, C, D     | Plot average |R| for each regulon source, for each cancer type and method |
| 5_2_make_dotplots.R                                 | Fig. 3A, B        | Plot alinearverage |R| for each cancer type using matched and mismatched regulons |
| misc_toy_plot_pos_vs_neg_corr.R                     | Fig. 4A           | Example plot illustrating type of correlation between activity and senstivity |
| 5_3_linear_modelling.R                              | Fig. S1           | Each term's contribution to the linear model (adjusted R-squared) |
| 6_1_make_faceted_barcherts_plots_signif_genes.R     | Fig. 4B, C, D, S3 | Faceted barchart for each cancer type with |R| stratified from 0.2 to 1 |
| 6_2_variance_per_gene_method_comparison.R           | Fig. S4           | Boxplot showing per gene variance across activity methods |
| 7_2_make_faceted_barcherts_cles.R                   | Fig. 5A, B, C, S5 | Faceted barchart for each cancer type with CLES (binary essentiality pred) |

