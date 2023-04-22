#Author: Emily Kibbler
#This code is to evaluate potential reference genes for follow-up studies by qPCR

#install.packages("readxl")
#install.packages("tidyverse")

library(readxl)
library(tidyverse)

all_data<-read_xlsx("GSE224615_DEGs.xlsx")
all_data<- all_data %>% rename("5019" = "5019-3") #strip the suffix from this column name
all_data<- all_data %>% rename("5057" = "5057-3") #same
all_data$AveExpr<-rowMeans(all_data[,5:40])

#filter data to include genes with log2fc <1 and average expression >10000 reads
ref_genes<-subset(all_data,all_data$log2FoldChange<1&all_data$log2FoldChange>-1&all_data$AveExpr>10000)
nrow(ref_genes) #454 candidate genes at this point
ref_genes$stdev<-NA
for(i in 1:nrow(ref_genes)){
  ref_genes$stdev[i]<-sd(ref_genes[i,5:40])
}
ref_genes$var<-ref_genes$stdev/ref_genes$AveExpr #standard deviation as a fraction of mean reads
#take the top 20 least variable genes
ref_genes<-ref_genes[order(ref_genes$var),][1:20,]

view(ref_genes[order(ref_genes$AveExpr,decreasing=FALSE),])


#Next steps: follow up in UCSC genome browser to further narrow down
#Start by investigating top candidate: HNRNPA3
#AKAP13: not great
#SET: seems like it would work too
#SLC44A2: another good candidate