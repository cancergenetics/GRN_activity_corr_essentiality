###This script processes aracne regulons

library(aracne.networks)
library(tidyverse)
library(binilib)
library(org.Hs.eg.db)

#prepare the regulons ----

regulons_list = c("regulonlaml", "regulongbm", "regulonluad", "reguloncoad", 
                  "regulonbrca", "regulonpaad", "regulonhnsc", "regulonstad",
                  "regulonblca", "regulonkirc")

regulon_cancer_aracne_dict <- c("regulonlaml" = "aml", "regulongbm" = "gbm", "regulonluad" =  "luad", 
                                "reguloncoad" = "coad", "regulonbrca" = "brca", "regulonpaad" = "paad", 
                                "regulonhnsc" = "hnsc", "regulonstad" = "stad", "regulonblca" = "blca", 
                                "regulonkirc" = "kirc")

for (regulon in regulons_list) {
  regulon.df <- as.data.frame(reg2tibble(eval(parse(text = regulon))), annotate = TRUE) # turn regulon object into df

  # ENTREZID -> HGNC_SYMBOL
  sources <- regulon.df$source
  sources_symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = sources, 
                                           columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")
  colnames(sources_symbols) <- c("source", "source_symbol")
  
  targets <- regulon.df$target      
  targets_symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = targets, 
                                           columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")
  colnames(targets_symbols) <- c("target", "target_symbol")
  
  regulon.df <- cbind(regulon.df, sources_symbols, targets_symbols) #merging
  
  #final network as a dataframe
  regulon.df <- regulon.df %>%
    dplyr::select("source_symbol", "target_symbol", "mor", "likelihood")
  colnames(regulon.df) <- c("tf", "target", "mor", "likelihood")
  
  #updated MOR to work for decoupleR 2.1
  regulon.df$mor <- sign(regulon.df$mor) * regulon.df$likelihood
  regulon.df <- regulon.df %>%
    subset(select = -likelihood) %>% #likelihood deprecated in decoupleR 2.1
    distinct(tf, target, .keep_all = TRUE) #duplicate edges not accepted
  
  #save regulons for further use
  write.csv(regulon.df, file = paste0("./regulons/aracne_regulons/", regulon_cancer_aracne_dict[regulon], "_aracne_regulon.csv"), row.names = FALSE)
  print(paste(regulon, "DONE"))
} 
