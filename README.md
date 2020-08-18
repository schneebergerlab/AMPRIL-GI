# AMPRIL-GI
The scripts used in the project "The evolutionary dynamics of genetic incompatibilities introduced by duplicated genes in Arabidopsis thaliana"

## genotyping

## GWAS
1) run gwas, e.g.: with the phenotype of presence of non-functionalized HPA copy 

Rscript gwas.LM.multipro.r gwas.data.Rdata HPA.phenotype.01.txt HPA 20

2) draw the mahattan plot with gwas output.

Rscript manhattan_plot.r HPA.gwas.result.txt HPA.gwas
