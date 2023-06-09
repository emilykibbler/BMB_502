#R code for analysis of: https://doi.org/10.1101/2023.02.09.527892
#Data is normalized reads of RNA seq on 23 patients with long COVID (LC) and 13 non-LC controls
#Author: Emily Kibbler
#The purpose of this code is to do a deep dive on the HLA genes
#install.packages("readxl)
#install.packages("tidyverse")

library(readxl)
library(tidyverse)

all_data<-read_xlsx("GSE224615_DEGs.xlsx")
all_data<- all_data %>% rename("5019" = "5019-3") #strip the suffix from this column name
all_data<- all_data %>% rename("5057" = "5057-3") #same

all_data$AveExpr<-rowMeans(all_data[,5:40])

#read in sample matrix
samples<-data.frame(read_xlsx("term_proj_sample_matrix.xlsx"))
#samples<-samples[,c("PID","lc_status","sex")]

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
#eliminate character string matches that are not HLA genes
HLA_genes<-subset(HLA_genes,HLA_genes$Gene.name!="HHLA3") 
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
  labs(fill="Avg read count")+
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

#change to boxplot view
#all the genes
HLA_genes %>% ggplot(aes(x=as.factor(Gene.name), y=Signal,color=lc_status)) + 
  geom_boxplot()+
  ylab("Normalized read count")+
  xlab("Gene")+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  ggtitle("HLA gene expression")

#all the genes, log2 scale
HLA_genes %>% ggplot(aes(x=as.factor(Gene.name), y=log2(Signal),color=lc_status)) + 
  geom_boxplot()+
  ylab("Normalized read count: log2 scale")+
  xlab("Gene")+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  ggtitle("HLA gene expression")

#Class II alone, log2 scale
subset(HLA_genes,HLA_genes$Class=="II") %>% ggplot(aes(x=as.factor(Gene.name), y=log2(Signal),color=lc_status)) + 
  geom_boxplot()+
  ylab("Normalized read count:log2 scale")+
  xlab("Gene")+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  labs(color="LC status")+
  ggtitle("HLA class II expression")
#ggsave("class_II_boxplot.png")


#view data table-wise for manual inspection

HLA_genes_table<-subset(all_data,grepl("HLA",all_data$Gene.name,fixed=FALSE)==TRUE)
HLA_genes_table<-subset(HLA_genes_table,HLA_genes_table$Gene.name!="HHLA3")
HLA_genes_table<-subset(HLA_genes_table,HLA_genes_table$Gene.name!="HHLA2")
HLA_genes_table$stdev<-NA
for(i in 1:nrow(HLA_genes_table)){
  HLA_genes_table$stdev[i]<-sd(HLA_genes_table[i,5:40])
}
HLA_genes_table$CV<-HLA_genes_table$stdev/HLA_genes_table$AveExpr

view(HLA_genes_table)
