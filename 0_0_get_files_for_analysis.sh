mkdir depmap_input_data
mkdir depmap_input_data/wrangled_input_data

mkdir tcga_data

mkdir regulons
mkdir ./regulons/aracne_regulons
mkdir ./regulons/grndb_regulons

mkdir results
mkdir ./results/results_matrices_per_combination
mkdir ./results/results_matrices_per_combination/aracne
mkdir ./results/results_matrices_per_combination/grndb
mkdir ./results/results_matrices_per_combination/dorothea

mkdir ./results/correlation_matrices
mkdir ./results/correlation_matrices/aracne
mkdir ./results/correlation_matrices/grndb
mkdir ./results/correlation_matrices/dorothea

mkdir ./results/barcharts_data
mkdir ./results/plots_data
mkdir ./results/go_enrichment

mkdir plots
mkdir ./plots/skittles_plots
mkdir ./plots/dotplots
mkdir ./plots/barcharts

cd ./tcga_data
wget https://xena.treehouse.gi.ucsc.edu/download/TumorCompendium_v10_PolyA_hugo_log2tpm_58581genes_2019-07-25.tsv
wget https://xena.treehouse.gi.ucsc.edu/download/clinical_TumorCompendium_v10_PolyA_2020-01-28.tsv
cd ..

#get depmap 21Q4 input data
cd ./depmap_input_data

#gene effect
wget -O CRISPR_gene_effect.csv https://ndownloader.figshare.com/files/31315996

#sample info
wget -O sample_info.csv https://ndownloader.figshare.com/files/31316011

#gene expression
wget -O CCLE_expression.csv https://ndownloader.figshare.com/files/31315882

cd ..

#grndb regulons
cd ./regulons/grndb_regulons
wget -O aml_TCGA-regulons.txt http://www.grndb.com/download/txt?condition=AML_TCGA
wget -O gbm_TCGA-regulons.txt http://www.grndb.com/download/txt?condition=GBM_TCGA
wget -O luad_TCGA-regulons.txt http://www.grndb.com/download/txt?condition=LUAD_TCGA
wget -O coad_TCGA-regulons.txt http://www.grndb.com/download/txt?condition=COAD_TCGA
wget -O brca_TCGA-regulons.txt http://www.grndb.com/download/txt?condition=BRCA_TCGA
wget -O paad_TCGA-regulons.txt http://www.grndb.com/download/txt?condition=PAAD_TCGA
wget -O hnsc_TCGA-regulons.txt http://www.grndb.com/download/txt?condition=HNSC_TCGA
wget -O stad_TCGA-regulons.txt http://www.grndb.com/download/txt?condition=STAD_TCGA
wget -O blca_TCGA-regulons.txt http://www.grndb.com/download/txt?condition=BLCA_TCGA
wget -O kirc_TCGA-regulons.txt http://www.grndb.com/download/txt?condition=KIRC_TCGA

cd ..
cd ..

#get master regulator list for go enrichment
wget https://www.cell.com/cms/10.1016/j.cell.2020.11.045/attachment/d212e665-9078-4ca1-b9c0-5f550380e158/mmc2.xlsx


#get decoupleR v2.1.8 from github for full code reproducibility
wget https://codeload.github.com/saezlab/decoupleR/zip/e8b80506c74e90d81bb51aa326fefd74422a3752
tar -czvf decoupleR.tar.gz e8b80506c74e90d81bb51aa326fefd74422a3752
