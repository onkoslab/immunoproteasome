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

#========Preparing dataframe for the genesets and extracting expression data for a set of genes ##########
dat = as.matrix(read_excel("data/r_input/immune_cells_gene_list.xlsx"))
x<-c(dat[,1])
geneSets <- list()
ttd<-list()

for(i in 1:18) {
  geneSets[[x[i]]] <-c(dat[i,4:ncol(dat)])# dat[,i]
  ttd<-append(ttd,c(dat[i,4:ncol(dat)]))
}

setwd("data/r_input/gene_exp")
aa<-list.files(path = "data/r_input/gene_exp",pattern = "csv$") 

for (i in 1:33) 
{
  cts <- as.matrix(read.csv(aa[i],row.names="Hybridization.REF"))
  include_list<-c(ttd)
  cts<-subset(cts, rownames(cts) %in% include_list)
  
  gsva_es <- gsva(cts, geneSets, mx.diff=1,kcdf="Poisson",method='gsva')
  tt<-strsplit(aa[i],split='_',fixed=TRUE)[[1]][3]
  tt<-strsplit(tt,split='.',fixed=TRUE)[[1]][1]
  
  ## writing GSVA score to output file
  setwd("data/r_output") 
  write.table(gsva_es, file=paste("GSVA_immune_cells_",tt,".tsv",sep =""), sep="\t")
 
 ### Designing and selecting the sample for different
  ttd1<-list()  
  setwd("data/r_input") 
  sample_id_info <- as.matrix(read.csv(paste('high_low_ip.csv',sep="")))
  sample_id_info_1<-sample_id_info[sample_id_info[,4]==tt, ]
  design <- cbind(sampleGroup1=1, sampleGroup2vs1=c(rep(0, nrow(sample_id_info_1[sample_id_info_1[,3]=='low', ])), rep(1, nrow(sample_id_info_1[sample_id_info_1[,3]=='high', ]))))
  ttd1<-append(ttd1,c(sample_id_info_1[,1])) #### get the sample id for the particular set of sample group
  
  gsva_es1<-subset(gsva_es, select=unlist(c(ttd1)))  ### get the gsva score for the particular set of sample groups
  tt<-strsplit(aa[i],split='_',fixed=TRUE)[[1]][3]
  tt<-strsplit(tt,split='.',fixed=TRUE)[[1]][1]
  
  
  #### Differential immune cells analysis
  fit <- lmFit(gsva_es1, design)
  fit <- eBayes(fit)
  
  setwd("data/r_output")
  write.csv(topTable(fit,number=Inf, coef="sampleGroup2vs1",resort.by="logFC"),file=paste('high_low_immuno_prtoeasome_Diff_immune_score_',tt,".csv",sep=""))
  setwd("data/r_input/gene_exp")
}