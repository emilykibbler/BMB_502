#R code for analysis of: https://doi.org/10.1101/2023.02.09.527892
#Data is normalized reads of RNA seq on 23 patients with long COVID (LC) and 13 non-LC controls
#Author: Emily Kibbler
#Purpose of this code:
#Generate pairwise scatter plots to see high level patterns of the data
#Generate a sample matrix report starting from the matrix in GEO and adding some information specific to our analysis


#install.packages("readxl)
#install.packages("tidyverse")

library(readxl)
library(tidyverse)

all_data<-read_xlsx("GSE224615_DEGs.xlsx")
all_data<- all_data %>% rename("5019" = "5019-3") #strip the suffix from this column name
all_data<- all_data %>% rename("5057" = "5057-3") #same

all_data$AveExpr<-rowMeans(all_data[,5:40]) #calculate average normalized reads for each gene in the data

#read in sample matrix
samples<-data.frame(read_xlsx("term_proj_sample_matrix.xlsx"))

#IDs to do pairwise scatterplots on
#three each of: female LC, female non-LC, male LC, male non-LC
small_list<-c("GSM7027483", "GSM7027491", "GSM7027494", "GSM7027479", "GSM7027508", "GSM7027484", "GSM7027501", "GSM7027486", "GSM7027503","GSM7027487", "GSM7027493", "GSM7027481") 
sample_subset<-samples%>%filter((Library.Name %in% small_list))
data_subset<-all_data[, which((names(all_data) %in% sample_subset$PID)==TRUE)]

#uncomment and run the following lines to see a large matrix with small panels
#pdf("pairwise.pdf")
#pairs(data_subset)
#dev.off()

#cut one sample from each condition to make each panel bigger
#two each: female LC, female non-LC, male LC, male non-LC
smaller_list<-c( "GSM7027491", "GSM7027494", "GSM7027508", "GSM7027484", "GSM7027486", "GSM7027503","GSM7027493", "GSM7027481") 
sample_subset2<-samples%>%filter((Library.Name %in% smaller_list))
data_subset2<-all_data[, which((names(all_data) %in% sample_subset2$PID)==TRUE)]

pdf("pairwise.pdf")
pairs(data_subset2)
dev.off()


#generate a sample matrix report to include in supplementary materials of our final paper
sample_subset2$pairwise_scatterplot<-"yes"
sample_subset$pairwise_scatterplot<-"no"
detailed_matrix<-rbind(sample_subset2,sample_subset[which(is.na(match(sample_subset$PID,sample_subset2$PID))),])
detailed_matrix$mapping<-"yes"
samples$pairwise_scatterplot<-"no"
samples$mapping<-"no"
detailed_matrix<-rbind(detailed_matrix,samples[which(is.na(match(samples$PID,detailed_matrix$PID))),])

#write.csv(detailed_matrix,"sample_matrix_table.csv",row.names=FALSE)
