#R code for analysis on https://doi.org/10.1101/2023.02.09.527892
#Work in progress
#Authors: Emily Kibbler and Catrina Spruce
#install.packages("readxl)
#install.packages("tidyverse")
#install.packages("ggrepel")

library(readxl)
library(tidyverse)
library(ggrepel)

all_data<-read_xlsx("GSE224615_DEGs.xlsx")
all_data<- all_data %>% rename("5019-3" = "5019") #strip the unnecessary suffix from this column name
all_data<- all_data %>% rename("5057-3" = "5057") #same

all_data$AveExpr<-rowMeans(all_data[,5:40])

#read in sample matrix
samples<-data.frame(read_xlsx("term_proj_sample_matrix.xlsx"))
#samples<-samples[,c("PID","lc_status","sex")]
view(samples[,c("PID","lc_status","sex")])

#IDs to do pairwise scatterplots on
small_list<-c("GSM7027483", "GSM7027491", "GSM7027494", "GSM7027479", "GSM7027508", "GSM7027484", "GSM7027501", "GSM7027486", "GSM7027503","GSM7027487", "GSM7027493", "GSM7027481") 
sample_subset<-samples%>%filter((Library.Name %in% small_list))
data_subset<-all_data[, which((names(all_data) %in% sample_subset$PID)==TRUE)]

#pdf("pairwise.pdf")
#pairs(data_subset)
#dev.off()

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


#MA plot
#migrated to separate code file

#manually look at most-changed genes
summary(all_data$log2FoldChange)
most_changed<-all_data[order(all_data$log2FoldChange),][1:10,] #ten most downregulated genes
most_changed<-rbind(most_changed,all_data[order(all_data$log2FoldChange, decreasing=TRUE),][1:10,])

view(most_changed)
view(all_data[order(all_data$padj),][1:20,])
view(subset(all_data,all_data$padj<0.05))

most_changed<-all_data[order(all_data$log2FoldChange),][1:250,] #most downregulated genes
most_changed<-rbind(most_changed,all_data[order(all_data$log2FoldChange, decreasing=TRUE),][1:250,])
#write.table(most_changed$ID,file="top500.txt",row.names=FALSE,col.names=FALSE)
#write.table(most_changed$ID,file="top500.txt",row.names=FALSE,col.names=FALSE)
#this can be used to feed into webgestalt or other

#candidates for best housekeeping genes
#pull from other code


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

