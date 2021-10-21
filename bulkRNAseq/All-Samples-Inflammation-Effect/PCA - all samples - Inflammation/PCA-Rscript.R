#------------------------------------- INFORMATION -------------------------------------------------

###RNAseq001 Felix###
#This transcriptomic data set is derived from mouse enriched adipocytes from male or female mice under 
#different stimulation conditions (non-inflamed (NI,water) vs inflammed (DSS)) and with/without Atg7 in adipocytes.
#This analysis is a priniciple component analysis (PCA) to cluster the samples and explain their variance
#Steps before: MultiQC on sequencing reads, kallisto, TXimport
#Script Date: October 2020 ###

#------------------------------------- DIRECTORY --------------------------------------------------  
#set working direcotry to folder where you want to execute your analysis
setwd("/Users/frichter/Desktop/RNAseqAd/Analysis")
getwd()

# If running for the first time you need to install these packages  
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.4")
BiocManager::install(version = "3.6")
BiocManager::install(version = "3.10")

BiocManager::install("tximport")
BiocManager::install("biomaRt")
BiocManager::install("DESeq2")
install.packages("ggplot2")

#----------------------------------- LOAD LIBRARIES -------------------------------------------- 

install.packages("ggplot2")
library(DESeq2)
library (ggplot2)
library(tximport)

# ----------------------------------- PRE DESEQ ----------------------------------------------------
#From the "TxImport" script you should have already created RData which has imported the kallisto data.
#During this previous step, the user is required to import the abundance.h5 files using tximport and saving the data (as RDS file) and values beforehand. 
#Using the RDS file, load in your txiKallisto object and metadatatable
txiKallistoObject <- readRDS("TxiKallistoObject_RNAseq_001")
metadata <- read.csv(file.choose())

# ----------------------------------- DESEQ: OVERALL OBJECT AFTER TXI----------------------------------------------------
# The uploaded RDS file should create a large list, which now needs to be converted into a DESeqDataSet.
# First, create an overall DESeq object from the txi.import containing all the data and variables present in the datafile.
# The argument "design" in this context are any variable that are present in the metadata depending on what you want to model, in my case probably all three major variables (condition, genotype, sex)
dds_object_all<-DESeqDataSetFromTximport(txiKallistoObject, colData=metadata, design=~condition+genotype+sex)

# ----------------------------------- PCA ----------------------------------------------------
# From Alina's script: Nick told her that the following code should also work: (counts per million: cpm <- log2cpm(gene_counts_10w)
# The rlog function allows me to transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size.
# For now I will do it as in Alina's original script and use the rlog function. 
normalised_all <- rlog(dds_object_all)

# Run prcomp with normalised overall data 
pcaData_all <- prcomp(t(assay(normalised_all)))
summary(pcaData_all)

# Next, we need to create PCA dataframe with PC1-PC4 and most important metadata for visualisation
pcaDataFrame_all <- data.frame(PC1=pcaData_all$x[,1],PC2=pcaData_all$x[,2],PC3=pcaData_all$x[,3],PC4=pcaData_all$x[,4],genotype=metadata$genotype,condition=metadata$condition,sex=metadata$sex)

# Install and load package: factoextra provides some easy-to-use functions to extract and visualize the output of multivariate data analyses
install.packages("factoextra")
library(factoextra)

# Install and load package: factoMineR is an exploratory data analysis methods to summarize, visualize and describe datasets (including PCA).
install.packages("FactoMineR")
library("FactoMineR")

# Scree plot of total dataset
# Eigenvalues correspond to the amount of the variation explained by each principal component (PC). This is part of the factoextra package.
# Explanation: You will get as many PC as samples in your dataset since only the same amount of PC as samples can explain all your variance. 
fviz_eig(pcaData_all)

# ----------------------------- Plot PCA -------------------------------------------------------
# Using ggplot function we can plot PCA data.

# Set colorvalues for plotting below
colorvalues <- c("NI" = "grey60", "DSS" = "steelblue")
colorvalues2 <- c("WT"= "grey60", "KO" = "steelblue")  

#ggplot of PC1+2, which according to the Scree plot shown above is explaining most of the variation in the data
ggplot(data=pcaDataFrame_all,aes(x=PC1,y=PC2,shape=sex, alpha=condition, col=genotype))+
  geom_point(aes(shape=sex, color=genotype), size=7)+scale_color_manual(values=colorvalues2)+
  labs(title="All samples",x="PC1(33.71%)", y="PC2(14.43%)", col="genotype")+
  theme(legend.title = element_text(color = "black", size = 15, face="bold"), legend.position = "bottom",legend.box = "vertical",legend.key.width = unit(1,"cm"), legend.text = element_text(color = "black", size = 15),  axis.title = element_text(size=15),  panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black"), legend.key = element_rect(colour = "white", fill = NA))+
  scale_shape_manual(values=c(15,17,8))+theme(plot.title = element_text(face="bold", size=25))+
  scale_x_continuous(limits = c(-50, 50))+
  scale_y_continuous(limits = c(-50,50))+
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))+
  geom_vline(xintercept=0, linetype="dotted", color="black", size=0.5)+
  geom_hline(yintercept=0, linetype="dotted", color="black", size=0.5)+
  scale_alpha_discrete(range=c(0.5,1))

ggplot(data=pcaDataFrame_all,aes(x=PC2,y=PC3,shape=sex, alpha=genotype, col=as.factor(condition)))+
  geom_point(aes(shape=sex, color=as.factor(condition)), size=7)+scale_color_manual(values=colorvalues)+
  labs(title="All samples",x="PC2(14.43%)", y="PC3(4.3%)", col="condition")+
  theme(legend.title = element_text(color = "black", size = 15, face="bold"), legend.position = "bottom",legend.box = "vertical",legend.key.width = unit(1,"cm"), legend.text = element_text(color = "black", size = 15),  axis.title = element_text(size=15),  panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black"), legend.key = element_rect(colour = "white", fill = NA))+
  scale_shape_manual(values=c(15,17,8))+theme(plot.title = element_text(face="bold", size=25))+
  scale_x_continuous(limits = c(-50, 50)) +scale_y_continuous(limits = c(-50, 50))+
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))+
  geom_vline(xintercept=0, linetype="dotted", color="black", size=0.5)+
  geom_hline(yintercept=0, linetype="dotted", color="black", size=0.5)+
  scale_alpha_discrete(range=c(0.5,1))
