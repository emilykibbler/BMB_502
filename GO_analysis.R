#install.packages("Bioconductor")
BiocManager::install("biomaRt")
library(biomaRt)


all_data<-read_xlsx("GSE224615_DEGs.xlsx")


most_changed<-all_data
most_changed$abs_change<-abs(most_changed$log2FoldChange)
most_changed<-most_changed[order(most_changed$abs_change,decreasing = TRUE),]
most_changed<-most_changed[1:250,] #subset on the top 250


ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
listFilters(ensembl)
mart=ensembl
#useMart(ensembl)

attributes_df<-as.data.frame(listAttributes(ensembl)) #what you want to retrieve
filters_df<-listFilters(ensembl) #how you want to find it

#find the correct keyword for what I want
view(subset(attributes_df,grepl("go",attributes_df$description,fixed=FALSE,ignore.case=TRUE)==TRUE))
view(subset(filters_df,grepl("ensembl",filters_df$name,fixed=FALSE,ignore.case=TRUE)==TRUE))


head(getBM(attributes=c('go_id','definition_1006','description','external_gene_name'),filters=	'ensembl_gene_id',most_changed$ID,mart=ensembl))

go_results<-getBM(attributes=c('go_id','definition_1006','description','external_gene_name'),filters=	'ensembl_gene_id',most_changed$ID,mart=ensembl)

