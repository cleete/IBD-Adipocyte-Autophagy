#Subset DEseq analysis for cohort/condition

##First perform the steps 1-2 to setup the right objects required for subsetting the cohorts.
## 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(ggplot2)

## ----------------------------------- DESEQ: SUBSETTING:CONDITION/COHORT ----------------------------------------

countdata = assay(dds_object_all)
gene_counts_NI <- countdata[ ,metadata$condition=="NI"]
gene_counts_DSS <- countdata[ ,metadata$condition=="DSS"]

coldata_NI = colData(dds_object_all) 
coldata_NI = coldata_NI[coldata_NI$condition=="NI",]

coldata_DSS = colData(dds_object_all) 
coldata_DSS = coldata_DSS[coldata_DSS$condition=="DSS",]

##----------------------------------- DESEQ: RUN: Non-inflamed cohort ----------------------------------------------------

#Create object
dds_object_NI= DESeqDataSetFromMatrix(countData =gene_counts_NI, colData = coldata_NI, design =~genotype+sex)

# Set level of comparison
dds_object_NI$sex<- relevel(dds_object_NI$condition, "M")
dds_object_NI$genotype<- relevel(dds_object_NI$genotype, "WT")

#Filter (optional)
keep <- rowSums(counts(dds_object_NI))>=10
dds_object_NI<- dds_object_NI[keep,]

#Run DESeq
dds_NI  = DESeq(dds_object_NI)

#list all coefficients
resultsNames(dds_NI)

#Save
save.image("Proj048.rda")

#Check if it worked 
resultsNames(dds_NI)

#Shrink results (choose one of the two)
results_NI<-lfcShrink(dds_NI,type="normal",coef="genotype_WT_vs_KO")
results_NI<-lfcShrink(dds_NI,type="apeglm",coef="genotype_KO_vs_WT")

#Order your results
results_dss_NI<- results_NI[order(results_NI$padj, decreasing = F),]

#Create dataframe and other columns
DESEQ_NI <- data.frame(results_NI)
DESEQ_NI$Genes <- row.names(DESEQ_NI) 
DESEQ_NI$Sig<-ifelse(DESEQ_NI$padj<0.05,"Yes","No")
DESEQ_NI$Labl<-ifelse(-log10(DESEQ_NI$padj)>4, DESEQ_NI$Genes, NA)

#Write as excel list
install.packages("xlsx")
library("xlsx")
write.csv(DESEQ_NI, "DEseq_NI_genotype.csv")

#Vulanco Plot
EnhancedVolcano(results_NI,
                lab = rownames(results_NI),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Non-inflamed: WT vs. KO',
                subtitle = NULL,
                xlim = c(-5,5),
                ylim= c(0,30),
                pCutoff = 5e-02,
                FCcutoff = 0,
                labSize = 4.0,
                labFace = "bold",
                colAlpha = 0.8,
                legendPosition = "none",
                col = c('grey', 'grey', "grey", "chartreuse4"))

#Plot values

plotCounts(dds_NI, gene = "Eda2r", intgroup = c("sex", "genotype"))


##----------------------------------- DESEQ: RUN: DSS cohort ----------------------------------------------------

#Create object
dds_object_DSS= DESeqDataSetFromMatrix(countData =gene_counts_DSS, colData = coldata_DSS, design =~genotype+sex)

# Set level of comparison
dds_object_DSS$sex<- relevel(dds_object_DSS$condition, "M")
dds_object_DSS$genotype<- relevel(dds_object_DSS$genotype, "WT")

#Filter (optional)
keep <- rowSums(counts(dds_object_DSS))>=10
dds_object_DSS<- dds_object_DSS[keep,]

#Run DESeq
dds_DSS = DESeq(dds_object_DSS)

#list all coefficients
resultsNames(dds_DSS)

#Save
save.image("Proj048.rda")

#Check if it worked 
resultsNames(dds_DSS)

#Shrink results (choose one of the two)
results_DSS<-lfcShrink(dds_DSS,type="normal",coef="genotype_WT_vs_KO")
results_DSS<-lfcShrink(dds_DSS,type="apeglm",coef="genotype_KO_vs_WT")

#Order your results
results_dss_DSS<- results_DSS[order(results_DSS$padj, decreasing = F),]

#Create dataframe and other columns
DESEQ_DSS <- data.frame(results_DSS)
DESEQ_DSS$Genes <- row.names(DESEQ_DSS) 
DESEQ_DSS$Sig<-ifelse(DESEQ_DSS$padj<0.05,"Yes","No")
DESEQ_DSS$Labl<-ifelse(-log10(DESEQ_DSS$padj)>4, DESEQ_DSS$Genes, NA)

#Write as excel list
install.packages("xlsx")
library("xlsx")
write.csv(DESEQ_DSS, "DEseq_DSS_genotype.csv")

#Vulanco Plot
EnhancedVolcano(results_DSS,
                lab = rownames(results_DSS),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'DSS: WT vs. KO',
                subtitle = NULL,
                xlim = c(-5,5),
                ylim= c(0,30),
                pCutoff = 5e-02,
                FCcutoff = 0,
                labSize = 4.0,
                labFace = "bold",
                colAlpha = 0.8,
                legendPosition = "none",
                col = c('grey', 'grey', "grey", "chartreuse4"))

#Plot values

plotCounts(dds_DSS, gene = "Adipoq", intgroup = c("genotype"))
