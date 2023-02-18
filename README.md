# Gene essentiality in cancer is better predicted by mRNA abundance than by gene regulatory network-inferred activity

The scripts available perform the analysis presented in the manuscript by Tudose, Bond and Ryan: Gene essentiality in cancer is better predicted by mRNA abundance than by gene regulatory network-inferred activity

Each script can be run in Rstudio. However, some of the scripts require more computational time, so it is recommended the analysis is at least partly run on a server using nohup and R CMD BATCH. This is especially the case for steps 1.0-4.0

IMPORTANT: Whilst most of the code is fully automated, the user should select at the start of each run the regulons for analysis (ARACNe, DoRothEA or GRNdb) or comment out the ones not necessary.

Scripts which need regulons selection:
*
*
*

Scripts which need data on all regulons in order to be able to run:
*
*
*

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
     [1] WebGestaltR_0.4.5      impute_1.68.0          org.Hs.eg.db_3.14.0    AnnotationDbi_1.56.2  
     [5] IRanges_2.28.0         S4Vectors_0.32.4       aracne.networks_1.20.0 viper_1.28.0          
     [9] Biobase_2.54.0         BiocGenerics_0.40.0    Cairo_1.6-0            cowplot_1.1.1         
    [13] reshape2_1.4.4         bmbstats_0.0.0.9001    decoupleR_2.1.8        ggrepel_0.9.3         
    [17] ggpubr_0.4.0           gplots_3.1.1           magrittr_2.0.2         forcats_1.0.0         
    [21] stringr_1.4.0          dplyr_1.1.0            purrr_1.0.1            readr_2.1.4           
    [25] tidyr_1.3.0            tibble_3.1.8           ggplot2_3.4.1          tidyverse_1.3.1       

    loaded via a namespace (and not attached):
      [1] colorspace_2.1-0       ggsignif_0.6.4         ellipsis_0.3.2         class_7.3-19          
      [5] XVector_0.34.0         fs_1.6.1               rstudioapi_0.14        proxy_0.4-27          
      [9] bit64_4.0.5            fansi_1.0.4            lubridate_1.9.2        xml2_1.3.3            
     [13] codetools_0.2-18       splines_4.1.0          doParallel_1.0.17      cachem_1.0.6          
     [17] jsonlite_1.8.4         apcluster_1.4.10       packrat_0.9.0          broom_1.0.3           
     [21] kernlab_0.9-32         dbplyr_2.3.0           png_0.1-8              compiler_4.1.0        
     [25] httr_1.4.4             backports_1.4.1        assertthat_0.2.1       Matrix_1.5-3          
     [29] fastmap_1.1.0          lazyeval_0.2.2         cli_3.6.0              htmltools_0.5.4       
     [33] tools_4.1.0            igraph_1.4.0           GenomeInfoDbData_1.2.7 gtable_0.3.1          
     [37] glue_1.6.2             doRNG_1.8.6            Rcpp_1.0.10            carData_3.0-5         
     [41] cellranger_1.1.0       Biostrings_2.62.0      vctrs_0.5.2            svglite_2.1.1         
     [45] nlme_3.1-152           iterators_1.0.14       rvest_1.0.3            timechange_0.2.0      
     [49] lifecycle_1.0.3        rngtools_1.5.2         gtools_3.9.4           rstatix_0.7.2         
     [53] zlibbioc_1.40.0        MASS_7.3-54            scales_1.2.1           hms_1.1.2             
     [57] parallel_4.1.0         memoise_2.0.1          segmented_1.6-2        stringi_1.7.12        
     [61] RSQLite_2.2.20         foreach_1.5.2          e1071_1.7-13           caTools_1.18.2        
     [65] GenomeInfoDb_1.30.1    systemfonts_1.0.4      rlang_1.0.6            pkgconfig_2.0.3       
     [69] bitops_1.0-7           lattice_0.20-44        htmlwidgets_1.6.1      bit_4.0.5             
     [73] tidyselect_1.2.0       plyr_1.8.8             R6_2.5.1               generics_0.1.3        
     [77] DBI_1.1.3              whisker_0.4.1          pillar_1.8.1           haven_2.5.1           
     [81] withr_2.5.0            mixtools_2.0.0         RCurl_1.98-1.10        survival_3.2-11       
     [85] KEGGREST_1.34.0        abind_1.4-5            modelr_0.1.10          crayon_1.5.2          
     [89] car_3.1-1              KernSmooth_2.23-20     utf8_1.2.3             plotly_4.10.1         
     [93] tzdb_0.3.0             grid_4.1.0             readxl_1.4.2           data.table_1.14.6     
     [97] blob_1.2.3             reprex_2.0.2           digest_0.6.31          munsell_0.5.0         
    [101] viridisLite_0.4.1     
