# GRN_activity_corr_essentiality

The scripts available perform the analysis presented in the manuscript by Tudose, Bond and Ryan: Gene essentiality in cancer is better predicted by mRNA abundance than by gene regulatory network-inferred activity

Each script can be run in Rstudio. However, some of the scripts require more computational time, so it is recommended the analysis is at least partly run on a server using nohup and R CMD BATCH. This is especially the case for steps 1.0-4.0

Whilst most of the code is fully automated, the user should select at the start of each run the regulons for analysis (ARACNe, DoRothEA or GRNdb) or comment out the ones not necessary.

For reproducibility, the package versions used in the original analysis are added via packrat. 
