#!/bin/bash
#this script creates ARACNe regulons from CCLE gene expression for BRCA and LUAD 
#NOTE: Aracne.jar can be downloaded from https://sourceforge.net/projects/aracne-ap/
java -Xmx24G -jar Aracne.jar -e ccle_regulons_input/luad_expr.txt -o aracne_output --tfs ccle_regulons_input/tfs_brca.txt --pvalue 1E-8 --seed 1 --calculateThreshold
java -Xmx24G -jar Aracne.jar -e ccle_regulons_input/luad_expr.txt -o aracne_output --tfs ccle_regulons_input/tfs_brca.txt --pvalue 1E-8 --seed 1
cd ../aracne_output
mv bootstrapNetwork_ul3atth75o35ngtur8ibskqq7s.txt aracne_luad_ccle.txt
cd ../ccle_regulons_input

java -Xmx24G -jar Aracne.jar -e ccle_regulons_input/brca_expr.txt -o aracne_output --tfs ccle_regulons_input/tfs_brca.txt --pvalue 1E-8 --seed 1 --calculateThreshold
java -Xmx24G -jar Aracne.jar -e ccle_regulons_input/brca_expr.txt -o aracne_output --tfs ccle_regulons_input/tfs_brca.txt --pvalue 1E-8 --seed 1
cd ../aracne_output
mv bootstrapNetwork_ul3atth75o35ngtur8ibskqq7s.txt aracne_brca_ccle.txt
cd ../ccle_regulons_input

java -Xmx24G -jar Aracne.jar -e ccle_regulons_input/coad_expr.txt -o aracne_output --tfs ccle_regulons_input/tfs_coad.txt --pvalue 1E-8 --seed 1 --calculateThreshold
java -Xmx24G -jar Aracne.jar -e ccle_regulons_input/coad_expr.txt -o aracne_output --tfs ccle_regulons_input/tfs_coad.txt --pvalue 1E-8 --seed 1
cd ../aracne_output
mv bootstrapNetwork_ul3atth75o35ngtur8ibskqq7s.txt aracne_coad_ccle.txt
cd ../ccle_regulons_input

java -Xmx24G -jar Aracne.jar -e ccle_regulons_input/paad_expr.txt -o aracne_output --tfs ccle_regulons_input/tfs_paad.txt --pvalue 1E-8 --seed 1 --calculateThreshold
java -Xmx24G -jar Aracne.jar -e ccle_regulons_input/paad_expr.txt -o aracne_output --tfs ccle_regulons_input/tfs_paad.txt --pvalue 1E-8 --seed 1
cd ../aracne_output
mv bootstrapNetwork_ul3atth75o35ngtur8ibskqq7s.txt aracne_paad_ccle.txt