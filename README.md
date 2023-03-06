# Immunoproteasome
Analysis of constitutive and Immunoproteasome expression level in solid tumors

Make the following directories in your current path to run the R scripts and notebooks
```
mkdir data
mkdir plots
mkdir supplementary_tables
mkdir data/tcga_tumor
mkdir data/gtex_normal
mkdir data/r_input
mkdir data/r_output
mkdir data/r_input/gene_exp
mkdir data/time_course
```
## Notebook details
### Run the notebooks and R scripts in the following order

  1. [Exploring_proteasome_expression_pan_cancer](https://github.com/Rahulncbs/Immunoproteasome/blob/main/Exploring_proteasome_expression_pan_cancer.ipynb)

 2. [immune_cells_gsva.R](https://github.com/Rahulncbs/Immunoproteasome/blob/main/immune_cells_gsva.R)

 3. [pathways_gsva.R](https://github.com/Rahulncbs/Immunoproteasome/blob/main/pathways_gsva.R)

 4. [epi_mes_gsva.R](https://github.com/Rahulncbs/Immunoproteasome/blob/main/epi_mes_gsva.R)

5. [Differential_immune_cells_enrichment](https://github.com/Rahulncbs/Immunoproteasome/blob/main/Differential_immune_cells_enrichment.ipynb)

6. [Differential_pathways_analysis.ipynb](https://github.com/Rahulncbs/Immunoproteasome/blob/main/Differential_pathways_analysis.ipynb)

7. [EMT_correlation_analysis.ipynb](https://github.com/Rahulncbs/Immunoproteasome/blob/main/EMT_correlation_analysis.ipynb)

8. [Copy_no_alteration_mutation_analysis.ipynb](https://github.com/Rahulncbs/Immunoproteasome/blob/main/Copy_no_alteration_mutation_analysis.ipynb)

9. [Survival_analysis_IP_CP.ipynb](https://github.com/Rahulncbs/Immunoproteasome/blob/main/Survival_analysis_IP_CP.ipynb)

10. [TNF_alpha_TGFB1_inducing_proteasome_expr_analysis.ipynb](https://github.com/Rahulncbs/Immunoproteasome/blob/main/TNF_alpha_TGFB1_inducing_proteasome_expr_analysis.ipynb)

11. [Methylation_and_upstream_pathways_IP_analysis.ipynb](https://github.com/Rahulncbs/Immunoproteasome/blob/main/Methylation_and_upstream_pathways_IP_analysis.ipynb)



**./data** contains all the input data for the analysis. tcga tumor gene expression data, GTEx normal gene expression data, tumor mutation burden, TCGA purity data,

**./plots** contains all the plots generated and provided in the figure

**./supplementart_tables** contains all the result

**./data/tcga_tumor** contains the TCGA gene expression data

**./data/gtex_normal** contains the GTEx gene expression data downloaded from https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/gtex_RSEM_Hugo_norm_count.gz 

**./data/r_input** contains the input data for computing gsva score

**./data/r_output** contains the resulting gsva scroe and differential enrichment score for immune cells and pathways

**./data/r_input/gene_exp** contains the gene expression from tcga to run the r script

**./data/time_course** contains the time course gene expression data which can be downloaded from here [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147405]





**Required data** 
```
All the data used in the analysis has been hosted at [link]. Please download the data, unzip it and put it in the data folder as mentioned below.
put all the data file from tcga_tumor folder to ./data/tcga_tumor directory.
put the gene signature file Epithelial_Mesenchymal_gene_list.csv, hall_mark_genes_df.csv, immune_and_pathways_gene_list.xlsx  to ./data/r_input directory.
leave the rest of the files in the ./data directory

Download the additonal GTEx normal gene expression data from  https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/gtex_RSEM_Hugo_norm_count.gz  and put it in the ./data/gtex_normal directory.
Download the time course data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147405  and put it in the ./data/time_course directory
```

### Required packages
#### Python 3.6.8
```
jupyter-notebook              5.7.11
matplotlib                    3.0.3
numpy                         1.16.3
seaborn                       0.9.0
scipy                         1.5.4
scanpy                        1.4.4
pandas                        0.25.0
statannot                    0.2.3
anndata                       0.6.22
lifelines                     0.26.3

```
#### R version 3.6.3
```
GSVA                          1.32.0
limma                         3.40.6
preprocessCore                1.46.0
GSEABase                      1.46.0
Biobase                       2.44.0
genefilter                    1.66.0
RColorBrewer                  1.1.3
readxl                        1.4.0
```
