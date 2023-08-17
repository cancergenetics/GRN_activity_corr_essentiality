#This script prepares the motif annotations to run SCENIC
#Follows the tutorial from SCENIC: http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html
#Output: processed motif annotation files for SCENIC

library(tidyverse)
library(SCENIC)
library(RcisTarget)
library(AUCell)

#motif annotations files
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather")

dir.create("cisTarget_databases"); setwd("cisTarget_databases")
for(featherURL in dbFiles) { #this may fail - in case that happens, manually download files from link above and save them in cisTarget_databases
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir 
}

setwd("..")
#modifying files metadata and saving them for use by SCENIC
db <- importRankings("./cisTarget_databases/hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather",
                     indexCol = "motifs")
names(db@rankings)[1] <- "features"
db@org <- "hgnc"
db@genome <- "hg19"
arrow::write_feather(db@rankings,
                     "./cisTarget_databases/hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather")

db <- importRankings("./cisTarget_databases/hg19-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather",
                     indexCol = "motifs")
names(db@rankings)[1] <- "features"
db@org <- "hgnc"
db@genome <- "hg19"
arrow::write_feather(db@rankings,
                     "./cisTarget_databases/hg19-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather")
