#R code for analysis on https://doi.org/10.1101/2023.02.09.527892
#Work in progress
#Authors: Emily Kibbler and Catrina Spruce
#install.packages("readxl)
#install.packages("tidyverse")
#install.packages("ggrepel")



library(readxl)
library(tidyverse)
library(ggrepel)
all_data<-read_xlsx("/Users/emilykibbler/Desktop/ekibbler/Documents/Personal/class/GSE224615_DEGs.xlsx")
all_data<- all_data %>% rename("5019-3" = "5019") #strip the suffix from this column name
all_data<- all_data %>% rename("5057-3" = "5057") #same

all_data$AveExpr<-rowMeans(all_data[,5:40])

#read in sample matrix
samples<-data.frame(read_xlsx("/Users/emilykibbler/Desktop/ekibbler/Documents/Personal/class/term_proj_sample_matrix.xlsx"))
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
HLA_genes<-subset(HLA_genes,HLA_genes$Gene.name!="HHLA3")
HLA_genes<-subset(HLA_genes,HLA_genes$Gene.name!="HHLA2")

#wondering if HLA gene expression changes a lot between people regardless of LC status or if it's tightly regulated in general

HLA_genes_table<-subset(all_data,grepl("HLA",all_data$Gene.name,fixed=FALSE)==TRUE)
HLA_genes_table<-subset(HLA_genes_table,HLA_genes_table$Gene.name!="HHLA3")
HLA_genes_table<-subset(HLA_genes_table,HLA_genes_table$Gene.name!="HHLA2")
HLA_genes_table$stdev<-NA
for(i in 1:nrow(HLA_genes_table)){
  HLA_genes_table$stdev[i]<-sd(HLA_genes_table[i,5:40])
}
HLA_genes_table$dev_over_mean<-HLA_genes_table$stdev/HLA_genes_table$AveExpr
view(HLA_genes_table)

#define class I and class II
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


HLA_genes%>% ggplot(aes(x=sex,y=Gene.name,fill=Signal))+
  geom_tile()+
  ylab("Gene name")+
  xlab("Sex")+
  scale_fill_gradient(low="white",high="purple")+
  ggtitle("HLA gene expression: Sex difference analysis")

subset(HLA_genes,HLA_genes$Class=="II") %>% ggplot(aes(x=lc_status,y=Gene.name,fill=Signal))+
  geom_tile()+
  ylab("Gene name")+
  xlab("LC status")+
  scale_fill_gradient(low="white",high="purple")+
  ggtitle("HLA class II gene expression in Long Covid")


subset(HLA_genes,HLA_genes$Class=="I") %>% ggplot(aes(x=as.factor(Gene.name), y=Signal)) + 
  geom_boxplot()+
  ylab("Normalized read count")+
  xlab("Gene")+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  ggtitle("HLA class I expression for all subjects")

subset(HLA_genes,HLA_genes$Class=="I"&HLA_genes$lc_status=="Non-LC") %>% ggplot(aes(x=as.factor(Gene.name), y=Signal)) + 
  geom_boxplot()+
  ylab("Normalized read count")+
  xlab("Gene")+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  ggtitle("HLA class I expression for non-LC subjects")

subset(HLA_genes,HLA_genes$Class=="II") %>% ggplot(aes(x=as.factor(Gene.name), y=log(Signal))) + 
  geom_boxplot()+
  ylab("Normalized reads: log10 scale")+
  xlab("Gene")+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  ggtitle("HLA class II expression for all subjects")

subset(HLA_genes,HLA_genes$Class=="II") %>% ggplot(aes(x=as.factor(Gene.name), y=Signal)) + 
  geom_boxplot()+
  ylab("Normalized reads")+
  xlab("Gene")+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  ggtitle("HLA class II expression for all subjects")

#MA plot

#all_data$AveExpr<-rowMeans(all_data[,5:40])

all_data$ave <- log2(rowMeans(all_data[, 5:40]))
ggplot(all_data, aes(x=ave, y=log2FoldChange, label=Gene.name)) +
  geom_point(aes(color = ifelse(pvalue<0.05, 'red', 'blue')),size=1) +
  geom_label_repel(data = dplyr::filter(all_data, abs(log2FoldChange) > 2),
                   size = 3,
                   box.padding = .5,
                   max.overlaps=25,
                   nudge_y = 1.0E-6) + ## This forces
  scale_colour_manual(labels = c("not sig", "p<0.05"), values=c('blue',    'red')) + 
  labs(color = "legend", 
       y = "log2FoldChange", 
       x = "Average log2 expression") +
  ggtitle("MA plot d.e. genes, LC vs non-LC") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12),
        panel.grid = element_line(linetype = "dashed", colour = "lightgrey"))

#Slightly modified plot with adjusted instead of raw p-values

all_data$ave <- log2(rowMeans(all_data[, 5:40]))
ggplot(all_data, aes(x=ave, y=log2FoldChange, label=Gene.name)) +
  geom_point(aes(color = ifelse(padj<0.05, 'red', 'blue')),size=1) +
  geom_label_repel(data = dplyr::filter(all_data, abs(log2FoldChange) > 2),
                   size = 3,
                   box.padding = .5,
                   max.overlaps=25,
                   nudge_y = 1.0E-6) + ## This forces
  scale_colour_manual(labels = c("not sig", "padj<0.05"), values=c('blue',    'red')) + 
  labs(color = "legend", 
       y = "log2FoldChange", 
       x = "Average log2 expression") +
  ggtitle("MA plot d.e. genes, LC vs non-LC") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12),
        panel.grid = element_line(linetype = "dashed", colour = "lightgrey"))

#manually look at most-changed genes
summary(all_data$log2FoldChange)
most_changed<-all_data[order(all_data$log2FoldChange),][1:10,] #ten most downregulated genes
most_changed<-rbind(most_changed,all_data[order(all_data$log2FoldChange, decreasing=TRUE),][1:10,])

view(most_changed)
view(all_data[order(all_data$padj),][1:20,])
view(subset(all_data,all_data$padj<0.05))

most_changed<-all_data[order(all_data$log2FoldChange),][1:250,] #most downregulated genes
most_changed<-rbind(most_changed,all_data[order(all_data$log2FoldChange, decreasing=TRUE),][1:250,])
write.table(most_changed$ID,file="top500.txt",row.names=FALSE,col.names=FALSE)
#this can be used to feed into webgestalt or other

#candidates for best housekeeping genes

ref_genes<-subset(all_data,all_data$pvalue>0.99&all_data$AveExpr>1000)
ref_genes$stdev<-NA
for(i in 1:nrow(ref_genes)){
  ref_genes$stdev[i]<-sd(ref_genes[i,5:40])
}
ref_genes<-subset(ref_genes,ref_genes$stdev<500)
view(ref_genes)
view(ref_genes$Gene.name)
