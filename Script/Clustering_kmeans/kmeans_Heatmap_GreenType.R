######################################################################
### kmeans clustering based on RNA-seq data
######################################################################

##Fig.5c, Fig.S2b

#library
library(ComplexHeatmap)
library(tidyverse)
library(viridis)
library(RColorBrewer)

#Current data
today <- Sys.Date()
print(today, quate=F)

###Data presprocessing----
#data input
#Gene expression data
count <- read.csv("ExpLevel/GeneExpressionMatrix.count.csv", row.names=1, header=T, sep="")
#Gene annotation Inf
GenesAno<-read.csv("GeneInfo/Lsativa_467_v5.annotation_info.txt", header=T, sep="\t")

#Genes with average counts <10 were excluded from the analysis
count <- count %>% 
  filter(rowMeans(count)>10)
#Arrangement geneID
str.split <- str_split(rownames(count), pattern = ".v5", simplify = TRUE)
geneID <- str_c(str.split[,1], str.split[,2], sep="_v5")
rownames(count) <- geneID

##TNN normalization
#Label for experiment
group <- factor(c(rep("01_White",42), rep("04_FL",42), rep("03_Blue",39), rep("02_Red",42)))

#Objectification of table
DGEList <- DGEList(counts = count, group = group)
#Design matrix of full model
design <- model.matrix(~ group)
#Calculation of Normalization (TMM coefficient) and Variance
DGEList <- calcNormFactors(DGEList, method="TMM") %>% 
  estimateGLMCommonDisp(design) %>% 
  estimateGLMTrendedDisp(design) %>% 
  estimateGLMTagwiseDisp(design)
#TMM normalized RPM 
TMMnorm <- cpm(DGEList) %>% as.data.frame()
TMMnorm.log <- log2(TMMnorm + 1)
#output TMMnorm
write.csv(TMMnorm.log, "ExpLevel/TMMnormalized_log2_count.csv", sep="")

###Statistic test by GLM (or exact test)----
#GLM
fit <- glmFit(DGEList, design)
#Comparing "group" (ANOVA-like)
lrt_gro <- glmLRT(fit, coef =c(2,3,4))
table_gro <- as.data.frame(topTags(lrt_gro, n = nrow(TMMnorm.log)))

sig <- table_gro %>% filter(FDR<0.001 | 
                              abs(logFC.group02_Red) >1 |
                              abs(logFC.group03_Blue) >1 |
                              abs(logFC.group04_FL) >1
)


TMMnorm.log.sig <- TMMnorm.log[rownames(sig),]

#MwtaData
MetaData <- read.csv("MetaData/SampleLabel.csv", header=T)
SampleNo <- paste("Sample_",formatC(1:56, width=2, flag="0"),"", sep="")

#Mean calculation in replicates
TMMnorm.Mean <- t(TMMnorm.log.sig) %>% 
  as.data.frame() %>%
  bind_cols(MetaData) %>% 
  pivot_longer(1:nrow(TMMnorm.log.sig), names_to="GeneID", values_to="ExpLevel") %>% 
  group_by(Light, Cultivar, Color, GeneID) %>% 
  summarize(ExpLevel=mean(ExpLevel)) %>% 
  pivot_wider(names_from=GeneID, values_from = ExpLevel) 

TMMnorm.Mean <- cbind.data.frame(SampleNo, TMMnorm.Mean)

#scaling
#Green types
TMMnorm.Mean <- TMMnorm.Mean %>% filter(Color=="Green")
TMMnorm.Mean.s <- scale(TMMnorm.Mean[5:ncol(TMMnorm.Mean)]) %>% 
  as.data.frame() 
rownames(TMMnorm.Mean.s) <-SampleNo[1:20]

##k-means
#Elbow method
#Compute and plot wss for k = 2 to k = 15.
#Green
#set parameter
n_clust <- 6 #Please decide this number by Gap or Elbow method
max_itr <- 100
k.max <- 20

wss <- sapply(1:k.max, 
              function(k){kmeans(t(TMMnorm.Mean.s), k, nstart = 50, iter.max = max_itr)$tot.withinss})

pdf(file = "kmeans/Fig/FigS2b_ElbowPlot_kmeans_Green.pdf", width=4.0, height=4.0)
plot(1:k.max, wss,las=3,
     type="b", pch = 1, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
dev.off()


#setseed
set.seed(29)
#k-means function
kmeans_out <- kmeans(t(TMMnorm.Mean.s), n_clust, nstart=50, iter.max=max_itr, algorithm = "Hartigan-Wong")
cl <- kmeans_out$cluster

#add cluster info to data matrix
TMMnorm.Mean.s.Clus <- TMMnorm.Mean.s %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(GeneID = colnames(TMMnorm.Mean.s)) %>% 
  mutate(clust = paste("Cluster_",formatC(cl, width=2,flag="0"), sep="")) 

#Matrix & label for heatmap  
Mat <- as.matrix(TMMnorm.Mean.s.Clus[1:nrow(TMMnorm.Mean.s.Clus),1:(ncol(TMMnorm.Mean.s.Clus)-2)])
#Light <- c(rep(1,14), rep(2,14), rep(3,14), rep(4,14))
Light <- c(rep(1,5), rep(2,5), rep(3,5), rep(4,5))

#Heatmap visualization by Complexheatmap
pdf(file ="kmeans/Fig/Fig5b_Heatmap_Exp_GreenType_k6.pdf",width=5, height=6)

Heatmap(Mat, 
        col=brewer.pal(11, "PRGn"),
        row_split = cl, 
        #row_km=10,
        #cluster_rows = FALSE, 
        #cluster_columns = FALSE,
        row_order=sort(colnames(TMMnorm.Mean.s)),
        #row_order=sort(names(cl)),
        column_order=sort(rownames(TMMnorm.Mean.s)),
        #column_split=4, 
        column_split=Light, 
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        show_row_names = FALSE,
        #column_dend_reorder = TRUE,
        top_annotation = HeatmapAnnotation(foo=anno_block(gp=gpar(fill=c("black", "red","blue","orange")),
                                                          labels=c("White","Red","Blue","FL"),
                                                          labels_gp=gpar(col="white", fontsize=10))),
        left_annotation = rowAnnotation(foo=anno_block(gp=gpar(fill=brewer.pal(6,"Dark2")),
                                                       labels=c(paste("GC",formatC(1:6, width=2, flag="0"),"", sep="")),
                                                       labels_gp=gpar(col="White", fontsize=10))),
        border = TRUE,
        use_raster = FALSE)

dev.off()

#Output of ExpTable in each cluster
GO <- read.csv("GeneInfo/GO_Info.csv",header=T)

for (i in formatC(1:6, width=2,flag="0")){
  print(i)
  ClusTable <- TMMnorm.Mean.s.Clus %>% filter(clust==paste("Cluster_",i,sep=""))
  ClusTable.GO <- left_join(ClusTable %>% select(clust,GeneID), GO, by="GeneID", na_matches = "never") %>% 
    na.omit()
  write.csv(ClusTable, paste("kmeans/ClusterInfo/Clus",i,"_Exp_GreenType.csv", sep=""))
  write.csv(ClusTable.GO, paste("kmeans/ClusterInfo/Clus",i,"_GO_GreenType.csv", sep=""), row.names = FALSE)
}


