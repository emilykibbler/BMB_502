ref_genes<-subset(all_data,all_data$pvalue>0.99&all_data$AveExpr>1500)
ref_genes$stdev<-NA
for(i in 1:nrow(ref_genes)){
  ref_genes$stdev[i]<-sd(ref_genes[i,5:40])
}

ref_genes<-subset(ref_genes,ref_genes$stdev<500)

ref_genes$CV<-NA
for(i in 1:nrow(ref_genes)){
  ref_genes$CV[i]<-(ref_genes[i,44]/ref_genes[i,42])
}

view(ref_genes)
view(ref_genes$Gene.name) 


#ref_gen_plot <- ggplot(ref_genes, aes(x = Gene.name, y = AveExpr, fill = 'deepskyblue4')) + 
#  geom_bar(position = position_dodge(width = 2.5), stat = "identity",fill = 'deepskyblue4' ) + 
# geom_errorbar(
#    aes(
#      x = Gene.name,
#      ymin = AveExpr - stdev,
#      ymax = AveExpr + stdev,
#    ),
#   width = 0.2,
#    position = position_dodge(width = 0.9),
#    stat = "identity"
#  ) +
#  ggtitle("Housekeeping Gene Candidates")
#print(ref_gen_plot)

#Comparison of common Ref Gene to LC associated genes
#common ref_genes
coref_genes<-subset(all_data,all_data$Gene.name == "GAPDH"|all_data$Gene.name == "ACTB")
coref_genes$stdev<-NA
for(i in 1:nrow(coref_genes)){
  coref_genes$stdev[i]<-sd(coref_genes[i,5:40])
}
coref_genes$CV<-NA
for(i in 1:nrow(coref_genes)){
  coref_genes$CV[i]<-(coref_genes[i,44]/coref_genes[i,42])
}

view(coref_genes)
view(coref_genes$Gene.name) 
#significant genes
sig_genes<-subset(all_data,all_data$Gene.name == "OR7D2"|all_data$Gene.name == "ALAS2")
sig_genes$stdev<-NA
for(i in 1:nrow(sig_genes)){
  sig_genes$stdev[i]<-sd(sig_genes[i,5:40])
}
sig_genes$CV<-NA
for(i in 1:nrow(sig_genes)){
  sig_genes$CV[i]<-(sig_genes[i,44]/sig_genes[i,42])
}

view(sig_genes)
view(sig_genes) 
#plot coref & sig_gene

coref_gen_plot <- ggplot(coref_genes, aes(x = Gene.name, y = as.numeric(CV), fill = 'deepskyblue4')) + 
  geom_bar(position = position_dodge(width = 0.4), stat = "identity",fill = 'deepskyblue4' ) + 
  labs(y= "Average Expr Var", x = "Gene")+
  ggtitle("Average Expression of Common Housekeeping Genes")
print(coref_gen_plot)

sig_gen_plot <- ggplot(sig_genes, aes(x = Gene.name, y = as.numeric(CV), fill = 'deepskyblue4')) + 
  geom_bar(position = position_dodge(width = 0.4), stat = "identity",fill = 'deepskyblue4' ) + 
  labs(y= "Average Expr Var", x = "Gene")+
  ggtitle("Long Covid Associated Gene")
print(sig_gen_plot)