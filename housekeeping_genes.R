#Author: Emily Kibbler
#This code is to evaluate potential housekeeping genes for follow-up studies by qPCR

#install.packages("readxl")
#install.packages("tidyverse")

library(readxl)
library(tidyverse)

all_data<-read_xlsx("GSE224615_DEGs.xlsx")
all_data<- all_data %>% rename("5019" = "5019-3") #strip the suffix from this column name
all_data<- all_data %>% rename("5057" = "5057-3") #same
all_data$AveExpr<-rowMeans(all_data[,5:40])

#filter data to include genes with log2fc <1 and average expression >10000 reads
HK_genes<-subset(all_data,all_data$log2FoldChange<1&all_data$log2FoldChange>-1&all_data$AveExpr>10000)
nrow(HK_genes) #454 candidate genes at this point
HK_genes$stdev<-NA
for(i in 1:nrow(HK_genes)){
  HK_genes$stdev[i]<-sd(HK_genes[i,5:40])
}
HK_genes$var<-HK_genes$stdev/HK_genes$AveExpr #standard deviation as a fraction of mean reads
#take the top 20 least variable genes
HK_genes<-HK_genes[order(HK_genes$var),][1:20,]

view(HK_genes[order(HK_genes$AveExpr,decreasing=FALSE),])


#Next steps: follow up in UCSC genome browser to further narrow down
#Start by investigating top candidate: HNRNPA3
#AKAP13: not great
#SET: seems like it would work too
#SLC44A2: another good candidate


#compare our selection methods to validated housekeeping genes
#source: /doi.org/10.1371/journal.pone.0260902
common_HKG<-c("ACTB", "B2M", "EF1α", "GAPDH", "GUSB", "HPRT", "PPIA", "RNA18S", "RPL13A", "TBP", "UBC", "YWHAZ")

documented_HKG<-all_data[which(!is.na(match(all_data$Gene.name,common_HKG))),]
documented_HKG$stdev<-NA
for(i in 1:nrow(documented_HKG)){
  documented_HKG$stdev[i]<-sd(documented_HKG[i,5:40])
}
documented_HKG$var<-documented_HKG$stdev/documented_HKG$AveExpr
view(documented_HKG)
