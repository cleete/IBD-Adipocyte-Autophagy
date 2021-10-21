#Data analysis including all samples to evaluate effect of colitis on adipocyte transcriptome
#The aim of this script is to analyze the dataset for differentially expressed genes based on homeostasis and DSS-treated group
#in visceral adipocytes. Thereby we regress other factors such as sex and genotype.

#------------------------------------ STEP1: Differential Gene Expression using DEseq2---------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install

BiocManager::install("apeglm")

install.packages("tximport")
library(tximport)
BiocManager::install("DESeq2")
library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggrepel)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

# ----------------------------------- STEP2: DESEQ: OVERALL OBJECT AFTER TXI----------------------------------------------------
# The uploaded RDS file should create a large list, which now needs to be converted into a DESeqDataSet.
# First, create an overall DESeq object from the txi.import containing all the data and variables present in the datafile.
# The argument "design" in this context are any variable that are present in the metadata depending on what you want to model, in my case probably all three major variables (condition, genotype, sex)


txiKallistoObject <- readRDS("TxiKallistoObject_RNAseq_001")
metadata <- read.csv(file.choose())

dds_object_all<-DESeqDataSetFromTximport(txiKallistoObject, colData=metadata, design=~condition+genotype+sex)

dds_object_all$condition <- factor(dds_object_all$condition, levels = c("NI","DSS"))


##----------------------------------- DESEQ: RUN: Water vs. DSS (looking at change that inflammation brings by) ----------------------------------------------------
## At this point we regress and genotype.


# Set level of comparison
dds_object_all$sex<- relevel(dds_object_all$sex, "M")
dds_object_all$condition<- relevel(dds_object_all$condition, "NI")

#Filter (optional after talking to Nick)
keep <- rowSums(counts(dds_object_all))>=10
dds_object_all<- dds_object_all[keep,]

#Run DESeq
dds_all = DESeq(dds_object_all)

#list all coefficients
resultsNames(dds_all)

#Save
save.image("Proj048_all replicates.rda")

#Shrink results (choose one of the two)
#results_all<-lfcShrink(dds_WT,type="normal",coef="condition_NI_vs_DSS")
results_all<-lfcShrink(dds_all,type="apeglm",coef="condition_DSS_vs_NI")

#Order your results
results_all<- results_all[order(results_all$padj, decreasing = F),]

#Create dataframe and other columns
DESEQ_all <- data.frame(results_all)
DESEQ_all$Genes <- row.names(DESEQ_all) 
DESEQ_all$Sig<-ifelse(DESEQ_all$padj<0.05,"Yes","No")
DESEQ_all$Labl<-ifelse(-log10(DESEQ_all$padj)>4, DESEQ_all$Genes, NA)

#Write as excel list
install.packages("xlsx")
library("xlsx")
write.csv(DESEQ_all, "211013_Final_DEseq_all_samples_Inflammation.csv")

#Vulanco Plot
EnhancedVolcano(results_all,
                lab = rownames(results_all),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'All Samples: NI vs. DSS',
                subtitle = NULL,
                xlim = c(-6,6),
                ylim= c(0,40),
                pCutoff = 5e-02,
                FCcutoff = 0,
                labSize = 4.0,
                labFace = "bold",
                colAlpha = 0.8,
                legendPosition = "none",
                col = c('grey', 'grey', "grey", "chartreuse4"))

#Plot values

plotCounts(dds_all, gene = "Rufy4", intgroup = c("condition"))


# ----------------------------------- Create normalized counts----------------------------------------------------

# Option 1: Get your normalized reads from DESeq objects
# Based on: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html 
dds <- estimateSizeFactors(dds_object_all)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
View(normalized_counts)
write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)
write.csv(normalized_counts, file="normalized_counts.csv")

