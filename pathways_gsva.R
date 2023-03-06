### Importing packages
library(limma)
library(GSVA)
library(preprocessCore)
library(GSEABase)
library(GSVAdata)
data(c2BroadSets) 
c2BroadSets
library(Biobase)
library(genefilter)
library(RColorBrewer)
library(readxl)

### Importing genes for 50 different hallmark pathways
dat = as.matrix(read.csv("/data/r_input/hall_mark_genes_df.csv", header = TRUE,row.names="X"))
geneSets <- list()
ttd<-list()
for(i in 1:50) 
{
  geneSets[[colnames(dat)[i]]] <- dat[,i]
  ttd<-append(ttd,dat[,i])
}

aa<-list.files(path = "data/r_input/gene_exp",pattern = "csv$") ## geting the gene exp file name
setwd("data/r_input/gene_exp")
for (i in 1:33)
{
  cts <- as.matrix(read.csv(aa[i],row.names="Hybridization.REF"))
  cts<-subset(cts, rownames(cts) %in% ttd)
  gsva_es <- gsva(cts, geneSets,kcdf="Poisson",mx.diff=1,method='gsva')
  tt<-strsplit(aa[i],split='_',fixed=TRUE)[[1]][3]
  tt<-strsplit(tt,split='.',fixed=TRUE)[[1]][1] 
  
  ### writing GSVA score to the output folder
  setwd("data/r_output")
  write.table(gsva_es, file=paste("GSVA_pathways_score_",tt,".tsv",sep =""), sep="\t")
  ######################################
  ttd1<-list()  
  setwd("data/r_input") 
  sample_id_info <- as.matrix(read.csv(paste('high_low_ip.csv',sep="")))
  sample_id_info_1<-sample_id_info[sample_id_info[,4]==tt, ]
  design <- cbind(sampleGroup1=1, sampleGroup2vs1=c(rep(0, nrow(sample_id_info_1[sample_id_info_1[,3]=='low', ])), rep(1, nrow(sample_id_info_1[sample_id_info_1[,3]=='high', ]))))
  ttd1<-append(ttd1,c(sample_id_info_1[,1])) #### get the sample id for the particular set of sample group
  
  gsva_es1<-subset(gsva_es, select=unlist(c(ttd1)))  ### get the gsva score for the particular set of sample groups
  
  
  #### Differential pathway analysis
  fit <- lmFit(gsva_es1, design) 
  fit <- eBayes(fit) 
  setwd("data/r_output")
  write.csv(topTable(fit,number=Inf, coef="sampleGroup2vs1",resort.by="logFC"),file=paste('high_low_immuno_prtoeasome_Diff_pathway_exp_',tt,".csv",sep=""))
  setwd("data/r_input/gene_exp") 
}
