## Load reauired packages

#install.packages("pheatmap")
library(tximport)
library(DESeq2)
library (BatchQC)
library("xlsx")
library(pheatmap)
library(dplyr)

# ----------------------------------- Create a DESeq object----------------------------------------------------

txiKallistoObject <- readRDS("TxiKallistoObject_RNAseq_001")
metadata <- read.csv(file.choose())
dds_object_all<-DESeqDataSetFromTximport(txiKallistoObject, colData=metadata, design=~condition+genotype+sex)
dds_object_all$genotype <- factor(dds_object_all$genotype, levels = c("WT","KO"))

# ----------------------------------- Create normalized counts----------------------------------------------------

# Option 1: Get your normalized reads from DESeq objects
# Based on: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html 
dds <- estimateSizeFactors(dds_object_all)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
View(normalized_counts)
write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)
write.csv(normalized_counts, file="normalized_counts.csv")

#Option 2: Make cpm file of your data (Not Working yet)
#countdata = assay(dds_object_all)
#counts.data_df <- as.data.frame(countdata)
#Gene_Counts_CPM <- log2CPM(counts.data_df)
#head(Gene_Counts_CPM)
#write.csv(Gene_Counts_CPM, file = "RNAseq_Counts_CPM") #save the CPM count



# -----------------------------------Read in normalised gene counts----------------------------------------------------

#Read in the normalised gene counts
Genes<- read.csv(file = '211013_Final_DEseq_all_samples_Inflammation_normCount.csv')

#Create dataframe 
Genes_df <- data.frame(Genes)

#Get rid of possible duplicates (tried and I have no duplicates)
Genes_df <-Genes_df[!duplicated(Genes_df$X), ]

#Get rid of N/As that might be problematic downstream 
Genes_df <-na.omit(Genes_df)

#Set my genenames as rownames 
rownames(Genes_df) <- Genes_df[,1]

#Delete column with the original gene names (1) as you do not need it anymore 
Genes_df <- Genes_df [, -1]


#Filter for the genes of interest
Genes_FFA_enrich <- as.data.frame(Genes_df [c("Acly","Cyp2f2","Acot3","Elovl4","Scd2","Cyp2e1","Acsm3","Fasn","Mid1ip1",
                                              "Acss2","Acaca","Aloxe3","Elovl6","Fads2","Pck1","Acot4","Fabp5","Ces1f","Etfbkmt",
                                              "Cd74","Ptges", "Lipg", "Cyp2d22","Irs2","Gk","Prkab2","Acot1","Ptgds","Acsm5","Akr1cl",
                                              "Pck2","Ptgs2","Hacd2","Acot11","Dgat2","Pdk2","Apoc2","Hsd17b10","Gpx1","Tnfrsf1a",
                                              "Insig1","Avpr1a","Sphk1","Cyp2u1","Ptgr1","Lpgat1","Ces1d","Pnpla3","Adipor2","Acss1",
                                              "Ephx2","Echdc2","H2-Ke6","Anxa1","Prkaa2","Gpat4","Cryl1","Pdk1","Nudt7","Acat2","Ptges2",
                                              "Adipoq","Appl2","Acsf2","Sirt1","Por","Prxl2b","Degs1","Elovl1","Adtrp","Echdc3","Pecr",
                                              "Dagla","Pam","Ltc4s","Fabp4","Dgat1","Acadsb","Cpt1a","Prkab1","Nucb2","Etfa","Eci1",
                                              "Eci1","Acox1","Ggt5","Mecr","Mtor","Slc27a4","Lipe","Abcd4","Slc25a17","Abcd3","Dld","Pibf1",
                                              "Sgpl1","Lpl","Lypla1","Qk","Cyb5a","Phyh","Eci3","Twist1","Hadha","Hsd17b4","Adipor1","Hadhb",
                                              "Pank2","Acsl5","Eci2","Akt2","Scp2","Ivd","Decr2","Lypla2","Etfdh","Nr1h2","Adh5","Gcdh",
                                              "Ech1","Prkag1","Auh"),])

Genes_Macroautophagy_enrich <- as.data.frame(Genes_df [c("Tlr2","Dcn","Rufy4","Bnip3","Retreg1","Kdr","Pik3c2b","Stbd1","Tigar",
                                                    "Lix1l","Lrrk2","Pim2", "Tbk1", "Ripk2", "Map1lc3b","Gabarapl1","Hmox1", "Atg14",
                                                    "Atp6v0a1","Tmem41b","Ctsd","Pink1","Nod1","Wdr45","Prkaa2","Ddrgk1","Bnip3l","Bag3",
                                                    "Trim13","Sirt1","Nprl3","Wdfy3","Sesn1","Mfsd8","Vps13c","Epg5","Htt","Prkn",
                                                    "Fez2","Ubxn2b","Huwe1","Uba5","Optn","Mtor","Htra2","Epm2a","Vps13d","Usp30",
                                                    "Hdac10","Atg16l1","Pikfyve","Ift88","Atg4b","Sh3glb1","Snx14","Map1lc3a","Pik3c2a",
                                                    "Ehmt2","Vamp8","Tomm7","Ubxn6","Wdr45b","Gba","Gabarap","Smurf1","Atp13a2","Atg13",
                                                    "Atg9a","Rab1b","Rab1a","Atg2a","Mtmr3","Tsc2","Tex264","Gnai3"),])



#Transpose your table to make matching it later for pheatmap easier 
Transposed_FFA_enrich <- as.data.frame(t(Genes_FFA_enrich))
Transposed_Autophagy_enrich <- as.data.frame(t(Genes_Macroautophagy_enrich))

#To properly use heatmap you should create a new dataframe (do not merge) that has the same row names as your filtered table. 
metadata_2 <- read.delim(file.choose())

#Set sample names as rownames 
rownames(metadata_2) <- metadata_2[,1]

#Create simpler annotations by choosing one column or now by subsetting to condition of interest. Pheatmap needs this annotation table.
annotation_df <-  as.data.frame(metadata[,9])
rownames(annotation_df) <- metadata[,1]

# Setup dataframe
Transposed_FFA_enrich <- as.data.frame(Transposed_FFA_enrich)
view(Transposed_FFA_enrich)

Transposed_Autophagy_enrich <- as.data.frame(Transposed_Autophagy_enrich)
view(Transposed_Autophagy_enrich)

#Get rid of Atg7Ad samples 
annotation_WT <- as.data.frame(metadata[-c(1,2,6,7,8,10,15,16,17,18,22,24),])

# Get rid of Atg7Ad samples in your subsetted gene tables 
Transposed_TNF_enrich <- as.data.frame(Transposed_TNF_enrich[-c(1,2,6,7,8,10,15,16,17,18,22,24),])
view(Transposed_TNF_enrich)

Transposed_Autophagy_enrich <- as.data.frame(Transposed_Autophagy_enrich[-c(1,2,6,7,8,10,15,16,17,18,22,24),])
view(Transposed_Autophagy_enrich)

FFA_enrich_num = as.matrix(Transposed_FFA_enrich[,1:122])
pheatmap(FFA_enrich_num)

Autophagy_enrich_num = as.matrix(Transposed_Autophagy_enrich[,1:75])
pheatmap(Autophagy_enrich_num)


FFA_enrich_num %>% 
  log2() -> FFA_enrich_num_log2
pheatmap(FFA_enrich_num_log2, scale="column")

Autophagy_enrich_num %>% 
  log2() -> Autophagy_enrich_num_log2
pheatmap(Autophagy_enrich_num_log2, scale="column")




FFA_enrich_logMeanSub <- FFA_enrich_num_log2 - colMeans((FFA_enrich_num_log2))
pheatmap(FFA_enrich_logMeanSub)

FFA_enrich_logMeanSub/colSds(as.matrix(FFA_enrich_num_log2)) ->
  FFA_enrich_zscore
pheatmap(TNF_enrich_zscore)

#To get to zscores, you can just take the log2 transformed dataset and use the scale argument to generate you z-score.

#Set your conditions as factor so that they appear in the order you want them to appear
annotation_WT$condition <- factor(annotation_WT$condition, levels=c("NI","DSS"))

#Create annotation colours list 
ann_colors_WT <- list(condition=c("NI"="gray90","DSS"="forestgreen"))

# Now you are ready to runn pheatmap. 
#Comment here: Genes are on the bottom and clustering on the side. This bugs me but I didnt have time to change yet. 
#7x9
Transposed_TNF_enrich %>% pheatmap()

pheatmap(FFA_enrich_num_log2, scale = "column", fontsize = 10)
pheatmap(Autophagy_enrich_num_log2, scale = "column", fontsize = 10)

#Order samples in a specific order
IdealOrder <- c("visceral.adipocyte.425", "visceral.adipocyte.426", "visceral.adipocyte.435", "visceral.adipocyte.432", "visceral.adipocyte.433", "visceral.adipocyte.441", "visceral.adipocyte.427", "visceral.adipocyte.428", "visceral.adipocyte.455", "visceral.adipocyte.430", "visceral.adipocyte.431", "visceral.adipocyte.439", "visceral.adipocyte.343", "visceral.adipocyte.344", "visceral.adipocyte.345", "visceral.adipocyte.350", "visceral.adipocyte.352", "visceral.adipocyte.359", "visceral.adipocyte.341", "visceral.adipocyte.347", "visceral.adipocyte.348", "visceral.adipocyte.342", "visceral.adipocyte.349", "visceral.adipocyte.351")
Autophagy_enrich_num_log2_ordered <- Autophagy_enrich_num_log2[IdealOrder,] 
pheatmap(Autophagy_enrich_num_log2_ordered, scale = "column", fontsize = 10,cluster_rows=F)

FFA_enrich_num_log2_ordered <- FFA_enrich_num_log2[IdealOrder,] 
pheatmap(FFA_enrich_num_log2_ordered, scale = "column", fontsize = 10,cluster_rows=F)
