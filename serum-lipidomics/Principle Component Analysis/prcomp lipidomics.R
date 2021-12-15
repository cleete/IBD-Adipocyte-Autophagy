## Lipidomics analysis using princomp or prcomp
# Both packages are built-in packages in R environment
install.packages("factoextra")
library(ggplot2)
library(factoextra)
library("FactoMineR")
library(tibble)
library(dplyr)

lipid.data <- read.csv("A3_pmol.csv", header=TRUE, row.names = 1)
lipid.per <- read.csv("A4_molper.csv", header=TRUE, row.names = 1)

#THIS IS AN IMPORTANT STEP: We are here removing all lipids which have a least NA (very strikt analysis)
lipid.data.new <- na.omit(lipid.data) # removing all NA values

#checking data quickly
head(lipid.data.new[, 1:6], 12)
summary(lipid.data.new)
view(lipid.data.new)

#Invert data so that samples are in rows and lipids are in columns
lipid.data.inv <- as.data.frame(t(lipid.data.new))

#Same thing can be done with the molpercentage table - Not done now
lipid.per <- na.omit(lipid.per)
head(lipid.per[, 1:6], 12)
summary(lipid.per)
lipid.per.inv <- as.data.frame(t(lipid.per))

#upload metadataset (needs to be same order as samples)
metadata <- read.csv("lipid_meta_inv.csv", header=TRUE, row.names = 1)

#combine lipid.data.new table and metadataset
lipid.data.new2 <- rbind(lipid.data.new, metadata[1:3,]) 

lipid.per2 <- rbind(lipid.per, metadata[1:3,])

#transpose dataset in a new dataframe
lipid.data.new3 = as.data.frame(t(lipid.data.new2))

lipid.per3 = as.data.frame(t(lipid.per2))

#Subsetting of lipid.data.new3 (e.g. by inflammatory condition)
lipid.ni <- lipid.data.new3[lipid.data.new3$condition == "NI",]
lipid.dss <- lipid.data.new3[lipid.data.new3$condition == "DSS",]

lipid.per.ni <- lipid.per3[lipid.per3$condition == "NI",]
lipid.per.dss <- lipid.per3[lipid.per3$condition == "DSS",]

#Use prcomp function to identify the different PCA components
PCA.all <- prcomp(lipid.data.inv, scale = TRUE, center = T)
summary(PCA.all)
#plot(my_pr, type = "l")
#biplot(my_pr, sclae = 0)

PCA.per.all <- prcomp(lipid.per.inv, scale= TRUE, center=T)
summary(PCA.per.all)

#PCA done on NI samples of pmol
lipid.ni_new <- lipid.ni %>% 
  mutate_at(vars(1:145),as.numeric)
PCA.ni <- prcomp(lipid.ni_new[,1:145], scale = TRUE, center = T)
summary(PCA.ni)

#PCA done only on NI samples of pmol%
lipid.per.ni_new <- lipid.per.ni %>% 
  mutate_at(vars(1:145),as.numeric)
PCA.per.ni <- prcomp(lipid.per.ni_new[,1:145], scale = TRUE, center = T)
summary(PCA.per.ni)

#PCA done on DSS samples of pmol
lipid.dss_new <- lipid.dss %>% 
  mutate_at(vars(1:145),as.numeric)
PCA.dss <- prcomp(lipid.dss_new[,1:145], scale = TRUE, center = T)
summary(PCA.dss)

#PCA done only on DSS samples of pmol%
lipid.per.dss_new <- lipid.per.dss %>% 
  mutate_at(vars(1:145),as.numeric)
PCA.per.dss <- prcomp(lipid.per.dss_new[,1:145], scale = TRUE, center = T)
summary(PCA.per.dss)


#Scree plot results
library(factoextra)
fviz_eig(PCA.all)
fviz_eig(PCA.ni)
fviz_eig(PCA.dss)

fviz_eig(PCA.per.all)
fviz_eig(PCA.per.ni)
fviz_eig(PCA.per.dss)


#Add PCA data to our dataset (pmol dataset)
lipid.pca.all <- cbind(lipid.data.new3, PCA.all$x[,1:2])
lipid.ni.pca <- cbind(lipid.ni, PCA.ni$x[,1:2])
lipid.dss.pca <- cbind(lipid.dss, PCA.dss$x[,1:2])

#Add PCA data to our dataset (pmol% dataset)
lipid.per.pca.all <- cbind(lipid.per3, PCA.per.all$x[,1:2])
lipid.per.ni.pca <- cbind(lipid.per.ni, PCA.per.ni$x[,1:2])
lipid.per.dss.pca <- cbind(lipid.per.dss, PCA.per.dss$x[,1:2])

# PLOT WITH GGPLOT
install.packages("ggplot2")
library(ggplot2)

colorvalues2 <- c("WT" = "grey60", "KO" = "steelblue")

#this PCA code works now really well and according to the colour schemes of my previous graphs
ggplot(lipid.pca.all, aes(x=PC1, y=PC2)) +
  scale_color_manual(values=colorvalues2) +
  geom_point(aes(shape=sex, color=condition), size=4, alpha=0.5)+
  labs(title="All samples",x="PC1(54.38%)", y="PC2(14.92%)")+
  theme(legend.title = element_text(color = "black", size = 15, face="bold"), text = element_text(size=20), legend.position = "bottom",legend.box = "vertical",legend.key.width = unit(1,"cm"), legend.text = element_text(color = "black", size = 15),  axis.title = element_text(size=15),  panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black"), legend.key = element_rect(colour = "white", fill = NA))+
  scale_shape_manual(values=c(15,17,8))+
  theme(plot.title = element_text(face="bold", size=25))+
  scale_x_continuous(limits = c(-35, 35))+
  scale_y_continuous(limits = c(-20,20))+
  geom_hline(yintercept=0,linetype='dashed')#+
  geom_vline(xintercept = 0, linetype='dashed')+
  scale_alpha_discrete(range=c(0.5,1))


ggplot(lipid.ni.pca, aes(x=PC1, y=PC2, col = genotype, fill = genotype))+
  geom_point(aes(shape=sex, color=genotype), size=4)+ 
  #stat_ellipse(geom = "polygon", level=0.95, alpha = 0.5)+
  labs(x="PC1(56.32%)", y="PC2(13.75%)")+
  scale_color_manual(values=colorvalues2) +
  theme(legend.title = element_text(color = "black", size = 15, face="bold"), legend.position = "bottom",legend.box = "vertical",legend.key.width = unit(1,"cm"), legend.text = element_text(color = "black", size = 15),  axis.title = element_text(size=15),  panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black"), legend.key = element_rect(colour = "white", fill = NA))+
  scale_shape_manual(values=c(15,17,8))+theme(plot.title = element_text(face="bold", size=25))+
  scale_x_continuous(limits = c(-40, 40))+
  scale_y_continuous(limits = c(-20,20))+
  geom_hline(yintercept=0,linetype='dashed')

# Ready to run: PCA for DSS colitis samples only
ggplot(lipid.dss.pca, aes(x=PC1, y=PC2))+
         scale_color_manual(values=colorvalues2) +
         geom_point(aes(shape=sex, color=genotype), size=4)+
         labs(title="DSS samples",x="PC1(65.77%)", y="PC2(10.24%)")+
         theme(legend.title = element_text(color = "black", size = 15, face="bold"), text = element_text(size=20), legend.position = "bottom",legend.box = "vertical",legend.key.width = unit(1,"cm"), legend.text = element_text(color = "black", size = 15),  axis.title = element_text(size=15),  panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black"), legend.key = element_rect(colour = "white", fill = NA))+
         scale_shape_manual(values=c(15,17,8))+
         theme(plot.title = element_text(face="bold", size=25))+
         scale_x_continuous(limits = c(-35, 35))+
         scale_y_continuous(limits = c(-20,20))+
         geom_hline(yintercept=0,linetype='dashed')+
         geom_vline(xintercept = 0, linetype='dashed') #+
        # scale_alpha_discrete(range=c(0.5,1))


## Plotting of pmol% results
ggplot(lipid.per.pca.all, aes(x=PC1, y=PC2, shape=condition)) +
  #scale_color_manual(values=colorvalues2) +
  geom_point(aes(shape=genotype, color=condition, alpha=sex), size=3, alpha=0.5)+
  #scale_color_manual(values=colorvalues2) +
  theme(legend.title = element_text(color = "black", size = 15, face="bold"), legend.position = "bottom",legend.box = "vertical",legend.key.width = unit(1,"cm"), legend.text = element_text(color = "black", size = 15),  axis.title = element_text(size=15),  panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black"), legend.key = element_rect(colour = "white", fill = NA))+
  scale_shape_manual(values=c(15,17,8))+theme(plot.title = element_text(face="bold", size=25))+
  scale_x_continuous(limits = c(-35, 35))+
  scale_y_continuous(limits = c(-20,20))+
  geom_hline(yintercept=0,linetype='dashed')
