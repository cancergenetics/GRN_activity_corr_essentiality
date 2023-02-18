install.packages("packrat")

packrat::init()

#install remotes to use the install_version function
install.packages("remotes")
library(remotes)

remotes::install_version("tidyverse", version = "1.3.1", repos = "http://cran.us.r-project.org")
remotes::install_version("magrittr", version = "2.0.2", repos = "http://cran.us.r-project.org") 
remotes::install_version("gplots", version = "3.1.1", repos = "http://cran.us.r-project.org")
remotes::install_version("ggpubr", version = "0.4.0", repos = "http://cran.us.r-project.org")
remotes::install_version("stringr", version = "1.4.0", repos = "http://cran.us.r-project.org") 
remotes::install_version("reshape2", version = "1.4.4", repos = "http://cran.us.r-project.org") 
remotes::install_version("cowplot", version = "1.1.1", repos = "http://cran.us.r-project.org") 
remotes::install_version("Cairo", version = "1.6-0", repos = "http://cran.us.r-project.org") 
remotes::install_version("ggrepel", version = "0.9.1", repos = "http://cran.us.r-project.org") 
remotes::install_version("WebGestaltR", version = "0.4.4", repos = "http://cran.us.r-project.org") 

#install bioconductor 3.14
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version = "3.14")

BiocManager::install("limma") #dependency for binilib
BiocManager::install("DESeq2") #dependency for binilib
BiocManager::install("fgsea") #dependency for binilib

remotes::install_github("fossbert/binilib")
remotes::install_github("mladenjovanovic/bmbstats")

BiocManager::install("aracne.networks")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("impute")

#install version 2.1.8 of decoupleR
install.packages("./decoupleR.tar.gz", repos = NULL, type ="source")
