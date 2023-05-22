############################################
## HCA based on transcriptome data
############################################

#Load library
library(tidyverse)
library(ggdendro)
library(edgeR)

#Current data
today <- Sys.Date()
print(today, quate=F)

#Input expression data matrix (RSEM_output)
count <- read.csv("GeneExp/Light_Cultivar/GeneExpressionMatrix.count.csv", header=T, sep="")

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

#Load SampleLabel information
MetaData <- read.csv("MetaData/SampleLabel.csv", header=T)

###HCA###
#Data scaling
TMMnorm.pca <- scale(t(TMMnorm.log))
#labeling
rownames(TMMnorm.pca)<-MetaData$label
#caluculate dist
TMMnorm.dist <- dist(TMMnorm.pca, method="euclidean")
#hclust
hca <- hclust(TMMnorm.dist, method="ward.D2")

#change format
dendro <- as.dendrogram(hca)
dendro.hca <- dendro_data(dendro, type = "rectangle")

#labeling
dendro.hca.labels <- dendro.hca$labels %>% left_join(MetaData, by="label")

#ggplot
ggplot(horiz = TRUE, segment(dendro.hca))+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dendro.hca.labels, aes(x, y, label = label, color = Light), size=3, hjust = 0, angle = 0)+
  ylim(-250, 800)+
  coord_flip()+ 
  scale_y_reverse(expand = c(0.2, 0))+
  theme(text=element_text(size=20),
        #axis.text.x = element_text(size=10, color="black", angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        )+
  scale_color_manual(values=c("azure4","red2","dodgerblue3","darkorange1"))+
  labs(y="Height")

ggsave("Fig/FigS8a_HCA_Color_Cultivar.pdf", width=15, height=20)


