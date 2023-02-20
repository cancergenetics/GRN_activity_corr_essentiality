# Gene essentiality in cancer is better predicted by mRNA abundance than by gene regulatory network-inferred activity

The scripts available perform the analysis presented in the manuscript by Tudose, Bond and Ryan: Gene essentiality in cancer is better predicted by mRNA abundance than by gene regulatory network-inferred activity

Each script can be run in Rstudio. However, some of the scripts require more computational time, so it is recommended the analysis is at least partly run on a server using nohup and R CMD BATCH. This is especially the case for steps 1.0-4.0

IMPORTANT: Whilst most of the code is fully automated, the user should select at the start of each run the regulons for analysis (ARACNe, DoRothEA or GRNdb) or comment out the ones not necessary.

Scripts which need regulons selection:
* 3_0_get_activity_values.R
* 4_0_run_correlations.R
* 5_0_make_skittlesplot.R
* 6_0_generate_data_faceted_barcharts_signif_drivers.R
* 6_1_make_faceted_barcharts_plots_signif_drivers.R
* 7_0_cles.R
* 7_1_generate_data_cles_barchart.R
* 7_2_make_faceted_barcharts_cles.R


Scripts which need data on all regulons in order to be able to run:
* 5_1_get_medians.R
* 5_3_linear_modelling.R

## For reproducibility
Step 1: Create a new Rproject in the directory where all the scripts are saved (R version 4.1.0, Rstudio version 1.4.1712)

Step 2: Run  0_1_create_packrat_session.R to install all packages  

For more information on packarat visit: https://github.com/rstudio/packrat/


## Overview figures from manuscript mapped to scripts

| Script                                              | Figures           | Brief description                                        |
|:----------------------------------------------------|:------------------|:---------------------------------------------------------|
| 2_1_make_barchart_cell_lines_per_cancer_type.R      | Fig. 2A           | Plot number of cell lines used for each cancer type in the CCLE |
| 5_0_make_skittlesplot.R                             | Fig. 2B, C, D     | Plot average \|R\| for each regulon source, for each cancer type and method |
| 5_2_make_dotplots.R                                 | Fig. 3A, B        | Plot alinearverage \|R\| for each cancer type using matched and mismatched regulons |
| misc_toy_plot_pos_vs_neg_corr.R                     | Fig. 4A           | Example plot illustrating type of correlation between activity and senstivity |
| 5_3_linear_modelling.R                              | Fig. S1           | Each term's contribution to the linear model (adjusted R-squared) |
| 6_1_make_faceted_barcherts_plots_signif_genes.R     | Fig. 4B, C, D, S3 | Faceted barchart for each cancer type with \|R\| stratified from 0.2 to 1 |
| 6_2_variance_per_gene_method_comparison.R           | Fig. S4           | Boxplot showing per gene variance across activity methods |
| 7_2_make_faceted_barcherts_cles.R                   | Fig. 5A, B, C, S5 | Faceted barchart for each cancer type with CLES (binary essentiality pred) |

## Session Info

    R version 4.1.0 (2021-05-18)
    Platform: x86_64-pc-linux-gnu (64-bit)
    Running under: Ubuntu 21.04

    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
    LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

    locale:
     [1] LC_CTYPE=en_IE.UTF-8       LC_NUMERIC=C               LC_TIME=en_IE.UTF-8       
     [4] LC_COLLATE=en_IE.UTF-8     LC_MONETARY=en_IE.UTF-8    LC_MESSAGES=en_IE.UTF-8   
     [7] LC_PAPER=en_IE.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
    [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_IE.UTF-8 LC_IDENTIFICATION=C       

    attached base packages:
    [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] WebGestaltR_0.4.4      impute_1.68.0          org.Hs.eg.db_3.14.0    AnnotationDbi_1.56.2  
     [5] IRanges_2.28.0         S4Vectors_0.32.4       aracne.networks_1.20.0 viper_1.28.0          
     [9] Biobase_2.54.0         BiocGenerics_0.40.0    Cairo_1.6-0            binilib_0.2.0         
    [13] cowplot_1.1.1          reshape2_1.4.4         bmbstats_0.0.0.9001    decoupleR_2.1.8       
    [17] ggrepel_0.9.1          ggpubr_0.4.0           gplots_3.1.1           magrittr_2.0.2        
    [21] forcats_1.0.0          stringr_1.4.0          dplyr_1.1.0            purrr_1.0.1           
    [25] readr_2.1.4            tidyr_1.3.0            tibble_3.1.8           ggplot2_3.4.1         
    [29] tidyverse_1.3.1        remotes_2.4.2         

    loaded via a namespace (and not attached):
      [1] readxl_1.4.2                backports_1.4.1             fastmatch_1.1-3            
      [4] systemfonts_1.0.4           igraph_1.4.0                plyr_1.8.8                 
      [7] lazyeval_0.2.2              splines_4.1.0               BiocParallel_1.28.3        
     [10] GenomeInfoDb_1.30.1         digest_0.6.31               foreach_1.5.2              
     [13] htmltools_0.5.4             fansi_1.0.4                 memoise_2.0.1              
     [16] doParallel_1.0.17           mixtools_2.0.0              tzdb_0.3.0                 
     [19] Biostrings_2.62.0           annotate_1.72.0             modelr_0.1.10              
     [22] matrixStats_0.63.0          svglite_2.1.1               timechange_0.2.0           
     [25] colorspace_2.1-0            blob_1.2.3                  rvest_1.0.3                
     [28] haven_2.5.1                 crayon_1.5.2                RCurl_1.98-1.10            
     [31] jsonlite_1.8.4              genefilter_1.76.0           iterators_1.0.14           
     [34] survival_3.2-11             glue_1.6.2                  gtable_0.3.1               
     [37] zlibbioc_1.40.0             XVector_0.34.0              DelayedArray_0.20.0        
     [40] car_3.1-1                   kernlab_0.9-32              apcluster_1.4.10           
     [43] abind_1.4-5                 scales_1.2.1                rngtools_1.5.2             
     [46] DBI_1.1.3                   rstatix_0.7.2               Rcpp_1.0.10                
     [49] viridisLite_0.4.1           xtable_1.8-4                bit_4.0.5                  
     [52] proxy_0.4-27                htmlwidgets_1.6.1           httr_1.4.4                 
     [55] fgsea_1.20.0                RColorBrewer_1.1-3          ellipsis_0.3.2             
     [58] pkgconfig_2.0.3             XML_3.99-0.13               dbplyr_2.3.0               
     [61] locfit_1.5-9.7              utf8_1.2.3                  tidyselect_1.2.0           
     [64] rlang_1.0.6                 munsell_0.5.0               cellranger_1.1.0           
     [67] tools_4.1.0                 cachem_1.0.6                cli_3.6.0                  
     [70] generics_0.1.3              RSQLite_2.3.0               broom_1.0.3                
     [73] fastmap_1.1.0               bit64_4.0.5                 fs_1.6.1                   
     [76] caTools_1.18.2              KEGGREST_1.34.0             packrat_0.9.0              
     [79] doRNG_1.8.6                 nlme_3.1-152                whisker_0.4.1              
     [82] xml2_1.3.3                  compiler_4.1.0              rstudioapi_0.14            
     [85] plotly_4.10.1               curl_5.0.0                  png_0.1-8                  
     [88] e1071_1.7-13                ggsignif_0.6.4              reprex_2.0.2               
     [91] geneplotter_1.72.0          stringi_1.7.12              lattice_0.20-44            
     [94] Matrix_1.5-3                vctrs_0.5.2                 pillar_1.8.1               
     [97] lifecycle_1.0.3             BiocManager_1.30.19         data.table_1.14.8          
    [100] bitops_1.0-7                GenomicRanges_1.46.1        R6_2.5.1                   
    [103] KernSmooth_2.23-20          gridExtra_2.3               codetools_0.2-18           
    [106] MASS_7.3-54                 gtools_3.9.4                assertthat_0.2.1           
    [109] SummarizedExperiment_1.24.0 DESeq2_1.34.0               withr_2.5.0                
    [112] GenomeInfoDbData_1.2.7      parallel_4.1.0              hms_1.1.2                  
    [115] grid_4.1.0                  class_7.3-19                segmented_1.6-2            
    [118] MatrixGenerics_1.6.0        carData_3.0-5               lubridate_1.9.2  
