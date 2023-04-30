#R code for analysis of: https://doi.org/10.1101/2023.02.09.527892
#Work in progress
#Author: Catrina Spruce, contributor: Emily Kibbler
#install.packages("readxl)
#install.packages("tidyverse")
#install.packages("ggrepel")


library(readxl)
library(tidyverse)
library(ggrepel)
## Load DEGs matrix from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224615
all_data<-read_xlsx("GSE224615_DEGs.xlsx")
all_data<- all_data %>% rename("5019" = "5019-3") #strip the suffix from this column name
all_data<- all_data %>% rename("5057" = "5057-3") #same
## Filter genes for expression >5.55 cpm.
all_data <- all_data %>% dplyr::filter(rowSums(all_data[, 5:40])>200)

## Check for normal read count distribution
hist(all_data$"5044", main="Histogram of expression, sample 5044",
     xlab="Expression level") 
## This is NOT normal! Huge tail
## TRy log2
hist(log2(all_data$"5044"), main="Histogram of log2 expression, sample 5044",
     xlab="Log2 expression level")
## Better. Let's check another
hist(log2(all_data$"5053"))
## Looks the same

## Log2 transform the dataset
## Separate out read counts:
counts <- all_data[ ,5:40]
counts <- log2(counts+1) ## Need to add +1 to avoid Inf values. log2(1) = 0
meta <- all_data %>% dplyr::select(Gene.name, ID, log2FoldChange, pvalue, padj)
all_data <- cbind(meta, counts)

## Load dataframe with sample metadata
samples<-data.frame(read_xlsx("C:/Users/cspruce/OneDrive - The Jackson Laboratory/Catrina/class/IntroBioinformatics/term_proj_sample_matrix.xlsx"))
samples <- samples %>% arrange(PID)

## Subset for 12 working samples
small_list<-c("GSM7027483", "GSM7027491", "GSM7027494", "GSM7027479", "GSM7027508", "GSM7027484", "GSM7027501", "GSM7027486", "GSM7027503","GSM7027487", "GSM7027493", "GSM7027481") 
sample_subset <- samples[(samples$Library.Name %in% small_list), ]
data_subset<-all_data[, which((colnames(all_data) %in% sample_subset$PID)==TRUE)]

## Perform PCA to look for separation of LC vs non-LC
#order <- order(colnames(all_data2)) ## Samples are already in numerical order
mat <- all_data[ , -c(1:5)]
pc <- prcomp( t(mat) ) # calculate PCA

# collect principle components 
counts.pca <- data.frame(cbind(pc$x, samples))
variance_exp <- (pc$sdev)^2 / sum(pc$sdev^2)

# plot
ggplot(data = counts.pca, aes(x = PC1, y = PC2, shape = as.factor(sex),  color=as.factor(lc_status))) +
  geom_point(size = 4) +
  scale_color_manual(breaks = c("LC", "Non-LC"), values=c("red", "green")) +
  labs(color = "Strain",
       shape = "Sex",
       x = paste0("PC1 variance explanded (", round(variance_exp[1]*100, 1), "%)"),
       y = paste0("PC2 variance explanded (", round(variance_exp[2]*100, 1), "%)"),
       title = "Clustering of samples by PCA") +
  labs(color = "Genotype",
       x = paste0("PC1 variance explanded (", round(variance_exp[1]*100, 1), "%)"),
       y = paste0("PC2 variance explanded (", round(variance_exp[2]*100, 1), "%)"),
       title = "Clustering of samples by PCA") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=10),
        panel.grid = element_line(linetype = "dashed", colour = "lightgrey"))

## There are two female samples that are way off. These are 5019 and 5057.
## Exclude poor samples and try again.
## Perform PCA to look for separation of LC vs non-LC
#order <- order(colnames(all_data2)) ## Samples are already in numerical order
mat2 <- all_data[ , -c(1:5)] %>% dplyr::select(-c("5019", "5057"))
samples2 <- samples %>% dplyr::filter(!PID %in% c("5019", "5057"))
pc <- prcomp( t(mat2) ) # calculate PCA

# collect principle components 
counts.pca <- data.frame(cbind(pc$x, samples2))
variance_exp <- (pc$sdev)^2 / sum(pc$sdev^2)

# plot
ggplot(data = counts.pca, aes(x = PC1, y = PC2, shape = as.factor(sex),  color=as.factor(lc_status))) +
  geom_point(size = 4) +
  scale_color_manual(breaks = c("LC", "Non-LC"), values=c("red", "green")) +
  labs(color = "Strain",
       shape = "Sex",
       x = paste0("PC1 variance explanded (", round(variance_exp[1]*100, 1), "%)"),
       y = paste0("PC2 variance explanded (", round(variance_exp[2]*100, 1), "%)"),
       title = "Clustering of samples by PCA") +
  labs(color = "Diagnosis",
       x = paste0("PC1 variance explanded (", round(variance_exp[1]*100, 1), "%)"),
       y = paste0("PC2 variance explanded (", round(variance_exp[2]*100, 1), "%)"),
       title = "Clustering of samples by PCA") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=10),
        panel.grid = element_line(linetype = "dashed", colour = "lightgrey"))

## This looks better! There is a general separation with LC on top and Non-LC on bottom.
## I'm surprised we don't see separation by sex, usually that is the first PC.
## It is tempting to cherry-pick the "best" LC and Non-LC samples and do differential analysis.
## There might be justification to do this if these LC samples had increased # of symptoms.
## And/or if the Non-LC didn't have a history of depression or other co-morbidities.
## We don't have the clinical data, so we can't.

## Which samples did we choose?? Color these purple
ggplot(data = counts.pca, aes(x = PC1, y = PC2, shape = as.factor(sex),  color=as.factor(lc_status))) +
  geom_point(size = 4) +
  scale_color_manual(breaks = c("LC", "Non-LC"), values=c("red", "green")) +
  geom_point(data=counts.pca[which((counts.pca$PID %in% sample_subset$PID)==TRUE), ], colour="purple", size=4) +
  labs(color = "Strain",
       shape = "Sex",
       x = paste0("PC1 variance explanded (", round(variance_exp[1]*100, 1), "%)"),
       y = paste0("PC2 variance explanded (", round(variance_exp[2]*100, 1), "%)"),
       title = "Clustering of samples by PCA") +
  labs(color = "Diagnosis",
       x = paste0("PC1 variance explanded (", round(variance_exp[1]*100, 1), "%)"),
       y = paste0("PC2 variance explanded (", round(variance_exp[2]*100, 1), "%)"),
       title = "Clustering of samples by PCA") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=10),
        panel.grid = element_line(linetype = "dashed", colour = "lightgrey"))

## So maybe better if we didn't pick outlier female on right.
## Could have picked 

## QC of p-values
hist(all_data$pvalue, main="Histogram of p-values",
     xlab="p-values")
## QC of adjusted p-values
hist(all_data$padj, main="Histogram of adjusted p-values",
     xlab="p-values")
## These look terrible. I would not trust any of these results

## MA plot of gene expression, with significantly d.e. genes colored red.
## Need to get average expression
all_data$ave <- rowMeans(all_data[, 6:41])
ggplot(all_data, aes(x=ave, y=log2FoldChange, label=Gene.name)) +
  geom_point(aes(color = ifelse(pvalue<0.05, 'red', 'blue')),size=1) +
  geom_label_repel(data = dplyr::filter(all_data, abs(log2FoldChange) > 1.5),
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

## What about coloring by adjusted p-value
ggplot(all_data, aes(x=ave, y=log2FoldChange, label=Gene.name)) +
  geom_point(aes(color = ifelse(padj<0.05, 'red', 'blue')),size=1) +
  geom_label_repel(data = dplyr::filter(all_data, padj<0.05),
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



## Try a volcano plot
ggplot(all_data, aes(x = log2FoldChange, y = -log(pvalue, 10), label = Gene.name)) +
    geom_point() +
    geom_point(aes(color = ifelse(pvalue<0.05, 'red', 'blue')),size=1) +
    geom_label_repel(data = dplyr::filter(all_data, abs(log2FoldChange)>1.5),
                     size = 3,
                     box.padding = .5,
                     max.overlaps=20) +
    labs(color = "legend", x = "log2FoldChange", y = "-log10(p-value)", color = "legend") +
    ggtitle("Volcano plot of d.e. genes, LC vs non-LC") +
  scale_colour_manual(labels = c("not sig", "p<0.05"), values=c('blue', 'red')) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12),
        panel.grid = element_line(linetype = "dashed", colour = "lightgrey"))

