#R code for analysis of: https://doi.org/10.1101/2023.02.09.527892
#Data is normalized reads of RNA seq on 23 patients with long COVID (LC) and 13 non-LC controls
#Author: Emily Kibbler
#This code is to evaluate potential housekeeping genes for follow-up studies by qPCR

#install.packages("readxl")
#install.packages("tidyverse")
library(readxl)
library(tidyverse)

all_data<-read_xlsx("GSE224615_DEGs.xlsx")
all_data$AveExpr<-rowMeans(all_data[,5:40])

#filter data to include genes with log2fc <1 and average expression >10000 reads
HK_genes<-subset(all_data,all_data$log2FoldChange<1&all_data$log2FoldChange>-1&all_data$AveExpr>10000)
nrow(HK_genes) #454 candidate genes at this point
HK_genes$stdev<-NA
for(i in 1:nrow(HK_genes)){
  HK_genes$stdev[i]<-sd(HK_genes[i,5:40])
}
HK_genes$CV<-HK_genes$stdev/HK_genes$AveExpr #standard deviation as a fraction of mean reads
#take the top 20 least variable genes
HK_genes<-HK_genes[order(HK_genes$CV),][1:20,]

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
documented_HKG$CV<-documented_HKG$stdev/documented_HKG$AveExpr
view(documented_HKG)


#candidates for best housekeeping genes
#pull from other code, analysis done by G. Ifijeh


#significantly expressed genes: plots
sig_genes<-subset(reformat,reformat$Gene.name=="OR7D2"|reformat$Gene.name=="ALAS2")
sig_genes$type<-"LC-associated"
doc<-subset(reformat,reformat$Gene.name=="GAPDH"|reformat$Gene.name=="ACTB")
doc$type<-"Documented HKG"
disc<-subset(reformat,reformat$Gene.name=="CBFA2T2"|reformat$Gene.name=="SP1")
disc$type<-"Candidate HKG"
df<-rbind(sig_genes,doc)
df<-rbind(df,disc)

a<-subset(df,df$lc_status=="Non-LC") %>% ggplot(aes(x=factor(Gene.name,level=c("GAPDH","ACTB","CBFA2T2","SP1","OR7D2","ALAS2")),y=log2(Signal),fill=type))+
  geom_boxplot()+
  ggtitle("Gene expression: Non-LC")+
  xlab("Gene name")+
  ylab("Normalized reads: log2 scale")+
  labs(fill="Gene type")+
  theme(axis.text.x=element_text(angle=45,hjust=1))
a

b<-subset(df,df$lc_status=="LC") %>% ggplot(aes(x=factor(Gene.name,level=c("GAPDH","ACTB","CBFA2T2","SP1","OR7D2","ALAS2")),y=log2(Signal),fill=type))+
  geom_boxplot()+
  ggtitle("Gene expression: LC")  +
  xlab("Gene name")+
  ylab("Normalized reads: log2 scale")+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  labs(fill="Gene type")

#install.packages("ggpubr")
library(ggpubr)

ggarrange(a,b,common.legend=TRUE)
#ggsave("genes.png")

c<-subset(df,df$lc_status=="Non-LC") %>% ggplot(aes(x=factor(Gene.name,level=c("GAPDH","ACTB","CBFA2T2","SP1","OR7D2","ALAS2")),y=log2(Signal),fill=type))+
  geom_violin()+
  ggtitle("Gene expression: Non-LC")+
  xlab("Gene name")+
  ylab("Normalized reads: log2 scale")+
  labs(fill="Gene type")+
  theme(axis.text.x=element_text(angle=45,hjust=1))


d<-subset(df,df$lc_status=="LC") %>% ggplot(aes(x=factor(Gene.name,level=c("GAPDH","ACTB","CBFA2T2","SP1","OR7D2","ALAS2")),y=log2(Signal),fill=type))+
  geom_violin()+
  ggtitle("Gene expression: LC")  +
  xlab("Gene name")+
  ylab("Normalized reads: log2 scale")+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  labs(fill="Gene type")

ggarrange(c,d,common.legend=TRUE)
#ggsave("violin.png")

