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

ref_genes<-subset(all_data,all_data$pvalue>0.99&all_data$AveExpr>1000)
ref_genes$stdev<-NA
for(i in 1:nrow(ref_genes)){
  ref_genes$stdev[i]<-sd(ref_genes[i,5:40])
}
ref_genes<-subset(ref_genes,ref_genes$stdev<500)

view(ref_genes)
view(ref_genes$Gene.name)

#Next steps: follow up in UCSC genome browser to further narrow down