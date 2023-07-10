############################################
## PCA based on RNA-seq
############################################

##Fig.5a

#library
library(tidyverse)
library(edgeR)

#Current data
today <- Sys.Date()
print(today, quate=F)

#Input expression data matrix (RSEM_output)
count <- read.csv("GeneExp/GeneExpressionMatrix.count.csv", header=T, sep="")

#Genes with average counts <10 were excluded from the analysis
count <- count %>% 
  filter(rowMeans(count)>10)

##TNN normalization
#Objectification of table & Calculation of TMM coefficient
DGEList <- DGEList(counts = count) %>% calcNormFactors(method="TMM")

#Pattern1----
#TMM normalized RPM 
TMMnorm <- cpm(DGEList) %>% as.data.frame()
TMMnorm.log <- log2(TMMnorm + 1)
#Arrangement geneID
str.split <- str_split(rownames(TMMnorm.log), pattern = ".v5", simplify = TRUE)
geneID <- str_c(str.split[,1], str.split[,2], sep="_v5")
rownames(TMMnorm.log) <- geneID

###PCA###
#Data scaling
TMMnorm.pca <- scale(t(TMMnorm.log))
#PCA by prcomp
pca <- prcomp(TMMnorm.pca)
summary(pca)

#Load SampleLabel information
Label <- read.csv("MetaData/SampleLabel.csv", header=T)
#Datamatrix arrangement for plot
pca.data <- as.data.frame(pca$x[, c(1:6)])
pca.data <- cbind.data.frame(pca.data, Label[,1:length(Label)])

#set turn
pca.data$Cultivar <- factor(pca.data$Cultivar, 
                         levels=c("01_RButt","02_RLeaf","03_ROak","04_RCos","05_FancR","06_HandR","07_RFire",
                                  "08_SunB","09_WineD","10_GOak","11_HandG","12_GLeaf","13_GButt","14_FrIce"))
                                  
pca.data$Light <- factor(pca.data$Light, 
                         levels=c("01_WhiteLED","02_RedLED","03_BlueLED","04_FL"))

#Propartion of variance in each pc
pc <- summary(pca)
pc1=format(signif(pc$importance[2]*100, digits = 3), nsmall=1) 
pc2=format(signif(pc$importance[5]*100, digits = 3), nsmall=1)
pc3=format(signif(pc$importance[8]*100, digits = 3), nsmall=1)
pc4=format(signif(pc$importance[11]*100, digits = 3), nsmall=1)
pc5=format(signif(pc$importance[14]*100, digits = 3), nsmall=1)
pc6=format(signif(pc$importance[18]*100, digits = 3), nsmall=1)


#Plot by ggplot2
#Light & Cultivar
ggplot(pca.data, aes(y=PC2, x=PC1, color=Light, shape=Cultivar))+
  guides(color=FALSE,fill=FALSE,shape=FALSE)+
  geom_point(size=3, alpha=0.9)+
  scale_shape_manual(values=c(0:9,15:18))+
  scale_color_manual(values=c("azure4","red2","dodgerblue3","darkorange1"))+
  #theme_cowplot(font_size = 16, line_size = 0.8)+
  theme_bw()+
  theme(text=element_text(size=14),axis.text=element_text(size=14,color="black"))+
  ylab(bquote("PC2 ("~.(pc2) ~"%)"))+
  xlab(bquote("PC1 ("~.(pc1) ~"%)"))

ggsave("Fig/Fig5a1_PCA_Light_Cultivar.pdf", width=3, height=3)


#Color & Cultivar
ggplot(pca.data, aes(y=PC2, x=PC1, color=Color, shape=Cultivar))+
  geom_point(size=3, alpha=0.9)+
  guides(color=FALSE,shape=FALSE)+
  scale_shape_manual(values=c(0:9,15:18))+
  scale_color_manual(values=c("green","deeppink"))+
  #theme_cowplot(font_size = 16, line_size = 0.8)+
  theme_bw()+
  theme(text=element_text(size=14),axis.text=element_text(size=14,color="black"))+
  ylab(bquote("PC2 ("~.(pc2) ~"%)"))+
  xlab(bquote("PC1 ("~.(pc1) ~"%)"))

ggsave("Fig/Fig5a2_PCA_Color_Cultivar.pdf", width=3, height=3)


