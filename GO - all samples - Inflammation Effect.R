#GO analysis according to Matthias' Script. For this I need a table with the log2-FC and all adj. p.values

#sometimes there can be issues loading particular packages, make sure that packages are in line with R versions.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

BiocManager::install("clusterProfiler", force = T)
library(clusterProfiler)
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

library("RColorBrewer")
library(ggplot2)
BiocManager::install("enrichplot", force = T)
library(enrichplot)


#Step1: load your csv file
deseq_all_inflammation <- read.csv("211013_Final_DEseq_all_samples_Inflammation.csv", header=T)

#Step2: curate the table
sign<-which(deseq_all_inflammation$padj<0.05) #defines that sign stand for all values which have a padj<0.05
deseq_all_inflammation_sign<-deseq_all_inflammation[sign,] #filter for only significant different genes (padj). This creates a shorter list of only significantly differentially expressed genes
deseq_all_inflammation_sign<-deseq_all_inflammation_sign[,c(1,3)] #since we don't need all columns, we can just get rid off all columns except gene name and FC

colnames(deseq_all_inflammation_sign)<-c("gene_name","log2FC") #We rename the columns
deseq_all_inflammation_sign$log2FC<-abs(deseq_all_inflammation_sign$log2FC) #For GO, we won't need to know in which direction a gene alters, hence we make all values absolute
deseq_all_inflammation_sign <- deseq_all_inflammation_sign[order(deseq_all_inflammation_sign$log2FC, decreasing=T),] #We order the values according to the FC
which(duplicated(deseq_all_inflammation_sign$gene_name)) #We confirm no dups in gene names

#Step 3: convert gene names to EntrezID
gene_symbols <- as.character(deseq_all_inflammation_sign[,1]) #We create a value containing only the names of the DEG
entrez_ids_GOI <- bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db") #convert gene symbols into EntrezID
##NOTE: 6.07% of input gene ID failed

which(duplicated(entrez_ids_GOI$ENTREZID)) #checking dups
deseq_all_inflammation_sign <- deseq_all_inflammation_sign[deseq_all_inflammation_sign$gene_name %in% entrez_ids_GOI$SYMBOL,]#keep only rows where entrezID was found
deseq_all_inflammation_sign[,1]<-entrez_ids_GOI$ENTREZID #assign entrez IDs to fold changes
colnames(deseq_all_inflammation_sign)<-c("entrez_id","log2FC")

#Step 4: prepare a list of background gene names
which(duplicated(deseq_all_inflammation$gene_name)) #confirming no dups in gene names
gene_symbols_bg<-as.character(deseq_all_inflammation[,1])
entrez_ids_bg<-bitr(gene_symbols_bg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db") #convert gene symbols into EntrezID
##NOTE: 15.59% of input genes failed to map

which(duplicated(entrez_ids_bg$ENTREZID)) #checking dups
deseq_all_inflammation <- deseq_all_inflammation[deseq_all_inflammation$X %in% entrez_ids_bg$SYMBOL,]#keep only rows where entrezID was found
deseq_all_inflammation[,1]<-entrez_ids_bg$ENTREZID #assign entrez IDs to fold changes
deseq_all_inflammation<-deseq_all_inflammation[,c(1,3)]
colnames(deseq_all_inflammation)<-c("entrez_id","log2FC")
deseq_all_inflammation$log2FC<-abs(deseq_all_inflammation$log2FC) #making FC absolute
deseq_all_inflammation <- deseq_all_inflammation[order(deseq_all_inflammation$log2FC, decreasing=T),]

#Step 4: GO annotation
goi_list<-deseq_all_inflammation_sign$entrez_id
bg_list<-deseq_all_inflammation$entrez_id
ego <- enrichGO(gene          = goi_list,
                universe      = bg_list,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
edox <- setReadable(ego, 'org.Mm.eg.db', 'ENTREZID')

#Check gene lists for the different GO terms
view(summary(ego))
summary(ego)

#Step 5:Visualization of GO analysis (Dot Plot)
coul <-colorRampPalette(brewer.pal(7,"RdYlBu"))(299)
dot <- dotplot(edox, showCategory=24, orderBy="p.adjust") + ggtitle("Water vs DSS")+ scale_color_gradientn(colors=coul) #chose own color gradient
coul2 <-colorRampPalette(rev(brewer.pal(7,"RdYlBu")))(299)
print(dot)

barplot(ego, showCategory=24)

#cnetplot, set number of pathways to include those of interest whilst maintaining legibility
install.packages("ggnewscale")
library(ggnewscale)

fc_vec<-deseq_all_inflammation_sign$log2FC
names(fc_vec)<-deseq_all_inflammation_sign$entrez_id
cnet <- cnetplot(edox, showCategory=10, foldChange=fc_vec)
print(cnet)
