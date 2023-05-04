#R code for analysis of: https://doi.org/10.1101/2023.02.09.527892
#Data is normalized reads of RNA seq on 23 patients with long COVID (LC) and 13 non-LC controls
#Author: Emily Kibbler
#Purpose of this code:
#Generate list(s) of top hits to run through String/DAVID/GO analysis using web tools

#install.packages("tidyverse")
#install.packages("readxl")
library(tidyverse)
library(readxl)

all_data<-read_xlsx("GSE224615_DEGs.xlsx")

most_changed<-all_data
most_changed$abs_change<-abs(most_changed$log2FoldChange)
most_changed<-most_changed[order(most_changed$abs_change,decreasing = TRUE),]
most_changed<-most_changed[1:250,] #subset on the top 250 most changed in either direction

#top 10 downregulated, top 10 upregulated
most_changed<-all_data
most_changed<-all_data[order(all_data$log2FoldChange),][1:10,] #most downregulated genes
most_changed<-rbind(most_changed,all_data[order(all_data$log2FoldChange, decreasing=TRUE),][1:10,])
most_changed<-most_changed[,c(1:4,41)]

view(most_changed)

#subset on fc more than 1
fc_more_than_one<-subset(most_changed,most_changed$abs_change>1)
#write.csv(fc_more_than_one,"fc_more_than_one.csv",row.names=FALSE)
