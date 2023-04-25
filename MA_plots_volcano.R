#Author: Catrina Spruce
#Contributor: Emily Kibbler

#install.packages("readxl")
#install.packages("tidyverse")
#install.packages("ggrepel")

library(readxl)
library(tidyverse)
library(ggrepel)
all_data<-read_xlsx("GSE224615_DEGs.xlsx")
all_data<- all_data %>% rename("5019" = "5019-3") #strip the suffix from this column name
all_data<- all_data %>% rename("5057" = "5057-3") #same

#plot using raw p values
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


#same plot, colors and labels adjusted to use adjusted p values
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


## Try a volcano plot
ggplot(all_data, aes(x = log2FoldChange, y = -log(pvalue, 10), label = Gene.name)) +
  geom_point() +
  geom_point(aes(color = ifelse(pvalue<0.05, 'red', 'blue')),size=1) +
  geom_label_repel(data = dplyr::filter(all_data, abs(log2FoldChange)>1.5),
                   size = 3,
                   box.padding = .5,
                   max.overlaps=20) +
  labs(x = "log2FoldChange", y = "-log10(p-value)", color = "legend") +
  ggtitle("Volcano plot of d.e. genes, LC vs non-LC")
theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
#
## Try a volcano plot with adj p values
ggplot(all_data, aes(x = log2FoldChange, y = -log(padj, 10), label = Gene.name)) +
  geom_point() +
  geom_point(aes(color = ifelse(padj<0.05, 'red', 'blue')),size=1) + +
  geom_label_repel(data = dplyr::filter(all_data, abs(log2FoldChange)>1.5),
                   size = 3,
                   box.padding = .5,
                   max.overlaps=20) +
  labs(x = "log2FoldChange", y = "-log10(p-value)", color = "legend") +
  ggtitle("Volcano plot of d.e. genes, LC vs non-LC")
theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
