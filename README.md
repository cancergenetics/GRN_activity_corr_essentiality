# Gene essentiality in cancer is better predicted by mRNA abundance than by gene regulatory network-inferred activity


<a href="https://zenodo.org/badge/latestdoi/599163601"><img src="https://zenodo.org/badge/599163601.svg" alt="DOI"></a>

The scripts available perform the analysis presented in the manuscript by Tudose, Bond and Ryan, available on *bioRxiv*: https://www.biorxiv.org/content/10.1101/2023.03.02.530664v2

Each script can be run in Rstudio. However, some of the scripts require more computational time, so it is recommended the analysis is at least partly run on a server using nohup and R CMD BATCH. This is especially the case for steps 1.0-4.0

IMPORTANT: Whilst most of the code is fully automated, the user should select at the start of each run the regulons for analysis (ARACNe, DoRothEA or GRNdb) or comment out the ones not necessary.

Scripts which need regulons selection:
* 3_0_get_activity_values.R
* 4_0_run_correlations.R
* 4_3_individual_genes_corr.R
* 5_0_make_skittlesplot.R
* 5_4_go_analysis.R
* 5_5_fishers_test_master_regulators.R
* 5_6_fishers_test_oncogenes.R
* 5_8_regulon_sizes_skittlesplots.R
* 5_9_unique_targets_skittlesplots.R
* 6_0_generate_data_faceted_barcharts_signif_drivers.R
* 6_1_make_faceted_barcharts_plots_signif_drivers.R
* 7_0_cles.R
* 7_1_generate_data_cles_barchart.R
* 7_2_make_faceted_barcharts_cles.R


Scripts which need data on all regulons in order to be able to run:
* 5_1_get_medians.R
* 5_3_linear_modelling.R
* 5_7_overlaps_strongly_correlated_and_GO.R
* 5_10_linear_model_regulon_size.R
* 5_11_linear_model_no_unique_targets.R
* 6_2_variance_per_gene_method_comparison.R 

## For reproducibility
Step 1: Create a new Rproject in the directory where all the scripts are saved (R version 4.1.0, Rstudio version 1.4.1712)

Step 2: Run  0_1_create_packrat_session.R to install all packages  

For more information on packarat visit: https://github.com/rstudio/packrat/


## Overview figures from manuscript mapped to scripts

| Script                                              | Figures           | Brief description                                        |
|:----------------------------------------------------|:------------------------|:---------------------------------------------------------|
| 2_1_make_barchart_cell_lines_per_cancer_type.R      | Fig. 2A                 | Plot number of cell lines used for each cancer type in the CCLE |
| 4_3_individual_genes_corr.R                         | Fig. S6B, C             | Plot correlations for individual genes with at least one R < -0.6 |
| 5_0_make_skittlesplot.R                             | Fig. 2B, C, D, S8A, S8B | Plot average \|R\| for each regulon source, for each cancer type and method |
| 5_2_make_dotplots.R                                 | Fig. 3A, B              | Plot alinearverage \|R\| for each cancer type using matched and mismatched regulons |
| misc_toy_plot_pos_vs_neg_corr.R                     | Fig. 4A                 | Example plot illustrating type of correlation between activity and senstivity |
| 5_3_linear_modelling.R                              | Fig. S1                 | Each term's contribution to the linear model (adjusted R-squared) |
| 5_7_overlaps_strongly_correlated_and_GO.R           | Fig. S5C, D, E, F, S6A  | Overlap strongly correlated genes and enriched GO terms overlaps |
| 5_8_regulon_sizes_skittlesplots.R                   | Fig. S2                 | Like 5_0, but stratified by regulon size |
| 5_9_unique_targets_skittlesplots.R                  | Fig. S3                 | Like 5_0, but stratified by number of unique targets |
| 5_10_linear_model_regulon_size.R                    | Fig. S4A                | Like 5_3, but including regulon size as a term |
| 5_11_linear_model_no_unique_targets.R               | Fig. S4A                | Like 5_3, but including number of unique targets as a term |
| 6_1_make_faceted_barcherts_plots_signif_genes.R     | Fig. 4B, C, D, S2       | Faceted barchart for each cancer type with \|R\| stratified from 0.2 to 1 |
| 6_2_variance_per_gene_method_comparison.R           | Fig. S3                 | Boxplot showing per gene variance across activity methods |
| 7_2_make_faceted_barcherts_cles.R                   | Fig. 5A, B, C, S4       | Faceted barchart for each cancer type with CLES (binary essentiality pred) |
