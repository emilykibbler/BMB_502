#R code for analysis on https://doi.org/10.1101/2023.02.09.527892
#Work in progress
#Author: Emily Kibbler
#install.packages("readxl)
#install.packages("tidyverse")
#install.packages("ggrepel")

library(readxl)
library(tidyverse)
library(ggrepel)

all_data<-read_xlsx("GSE224615_DEGs.xlsx")
all_data<- all_data %>% rename("5019-3" = "5019") #strip the suffix from this column name
all_data<- all_data %>% rename("5057-3" = "5057") #same

all_data$AveExpr<-rowMeans(all_data[,5:40])

#read in sample matrix
samples<-data.frame(read_xlsx("term_proj_sample_matrix.xlsx"))
#samples<-samples[,c("PID","lc_status","sex")]

#IDs to do pairwise scatterplots on
small_list<-c("GSM7027483", "GSM7027491", "GSM7027494", "GSM7027479", "GSM7027508", "GSM7027484", "GSM7027501", "GSM7027486", "GSM7027503","GSM7027487", "GSM7027493", "GSM7027481") 
sample_subset<-samples%>%filter((Library.Name %in% small_list))
data_subset<-all_data[, which((names(all_data) %in% sample_subset$PID)==TRUE)]

pdf("pairwise.pdf")
pairs(data_subset)
dev.off()

#head(all_data)
reformat<-data.frame(matrix(ncol=6,nrow=0))
colnames(reformat)<-c("PID","log2FoldChange","pvalue","padj","Signal","Gene.name")

#put data in single column to make it more ggplot friendly
for(i in 1:36){
  temp<-all_data[,c(1:4,(4+i),41)]
  temp$PID<-colnames(temp)[5]
  colnames(temp)[5]<-"Signal"
  temp<-merge(temp,samples,by="PID")
  reformat<-rbind(reformat,temp)
}


HLA_genes<-subset(reformat,grepl("HLA",reformat$Gene.name,fixed=FALSE)==TRUE)
HLA_genes<-subset(HLA_genes,HLA_genes$Gene.name!="HHLA3") #character string matches that are not HLA genes
HLA_genes<-subset(HLA_genes,HLA_genes$Gene.name!="HHLA2")



#label class I and class II
HLA_genes$Class<-NA
for(i in 1:nrow(HLA_genes)){
  if(HLA_genes$Gene.name[i]=="HLA-A"|HLA_genes$Gene.name[i]=="HLA-B"|HLA_genes$Gene.name[i]=="HLA-C"){
    HLA_genes$Class[i]<-"I"}else{HLA_genes$Class[i]<-"II"}
}


HLA_genes%>% ggplot(aes(x=lc_status,y=Gene.name,fill=Signal))+
  geom_tile()+
  ylab("Gene name")+
  xlab("LC status")+
  scale_fill_gradient(low="white",high="purple")+
  ggtitle("HLA gene expression: Long Covid analysis")
#ggsave("HLA_heatmap.png")



HLA_genes%>% ggplot(aes(x=sex,y=Gene.name,fill=Signal))+
  geom_tile()+
  ylab("Gene name")+
  xlab("Sex")+
  scale_fill_gradient(low="white",high="purple")+
  ggtitle("HLA gene expression: Sex difference analysis")
#unremarkable differences between sexes

#looking at class II alone (the lower-expressed class) may make changes more visible
subset(HLA_genes,HLA_genes$Class=="II") %>% ggplot(aes(x=lc_status,y=Gene.name,fill=Signal))+
  geom_tile()+
  ylab("Gene name")+
  xlab("LC status")+
  scale_fill_gradient(low="white",high="purple")+
  ggtitle("HLA class II gene expression in Long Covid")
#ggsave("HLA_classII_heatmap.png")

#look at class I alone to be fair
subset(HLA_genes,HLA_genes$Class=="I") %>% ggplot(aes(x=lc_status,y=Gene.name,fill=Signal))+
  geom_tile()+
  ylab("Gene name")+
  xlab("LC status")+
  scale_fill_gradient(low="white",high="purple")+
  ggtitle("HLA class I gene expression in Long Covid")
#unremarkable


#Investigate whether HLA genes are tightly regulated
#i.e., not differential expressed between people, or are differentially expressed but not due to LC

HLA_genes %>% ggplot(aes(x=as.factor(Gene.name), y=Signal)) + 
  geom_boxplot()+
  ylab("Normalized read count")+
  xlab("Gene")+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  ggtitle("HLA expression for all subjects")
#ggsave("HLA_expression_summary_boxplot.png")

HLA_genes %>% ggplot(aes(x=as.factor(Gene.name), y=Signal,color=lc_status)) + 
  geom_boxplot()+
  ylab("Normalized read count")+
  xlab("Gene")+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  ggtitle("HLA gene expression")
#ggsave("HLA_gene_expression_boxplot.png")

#HLA genes for everyone, log scale

HLA_genes %>% ggplot(aes(x=as.factor(Gene.name), y=log2(Signal))) + 
  geom_boxplot()+
  ylab("Normalized read count:log2 scale")+
  xlab("Gene")+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  ggtitle("HLA expression for all subjects")


#class I and class II are different expression levels
#split into separate graphs to see better
subset(HLA_genes,HLA_genes$Class=="I") %>% ggplot(aes(x=as.factor(Gene.name), y=Signal)) + 
  geom_boxplot()+
  ylab("Normalized read count")+
  xlab("Gene")+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  ggtitle("HLA class I expression for all subjects")

subset(HLA_genes,HLA_genes$Class=="II"&HLA_genes$lc_status=="Non-LC") %>% ggplot(aes(x=as.factor(Gene.name), y=Signal)) + 
  geom_boxplot()+
  ylab("Normalized read count")+
  xlab("Gene")+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  ggtitle("HLA class II expression for non-LC subjects")

subset(HLA_genes,HLA_genes$Class=="II") %>% ggplot(aes(x=as.factor(Gene.name), y=Signal)) + 
  geom_boxplot()+
  ylab("Normalized reads")+
  xlab("Gene")+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  ggtitle("HLA class II expression for all subjects")

#put Y axis in log scale
subset(HLA_genes,HLA_genes$Class=="II") %>% ggplot(aes(x=as.factor(Gene.name), y=log(Signal))) + 
  geom_boxplot()+
  ylab("Normalized reads: log10 scale")+
  xlab("Gene")+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  ggtitle("HLA class II expression for all subjects")





#view data table-wise

HLA_genes_table<-subset(all_data,grepl("HLA",all_data$Gene.name,fixed=FALSE)==TRUE)
HLA_genes_table<-subset(HLA_genes_table,HLA_genes_table$Gene.name!="HHLA3")
HLA_genes_table<-subset(HLA_genes_table,HLA_genes_table$Gene.name!="HHLA2")
HLA_genes_table$stdev<-NA
for(i in 1:nrow(HLA_genes_table)){
  HLA_genes_table$stdev[i]<-sd(HLA_genes_table[i,5:40])
}
HLA_genes_table$var<-HLA_genes_table$stdev/HLA_genes_table$AveExpr

view(HLA_genes_table)