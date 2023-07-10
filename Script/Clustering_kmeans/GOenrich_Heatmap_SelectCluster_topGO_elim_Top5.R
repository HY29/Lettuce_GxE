######################################################################
### Heatmap visualization of GO enrichment results
######################################################################

##Fig.5e

#library
library(tidyverse)
library(ComplexHeatmap)

#Current data
today <- Sys.Date()
print(today, quate=F)

#Data input
#Red
Clus01_Red <- read.csv("kmeans/ClusterInfo/topGO/Clus01_RedType_GenTable_classic_FDR005.csv",header=T) %>% 
  select(GO.ID) %>% 
  slice(1:3) %>% 
  mutate("Clus"="RC01")
Clus02_Red <- read.csv("kmeans/ClusterInfo/topGO/Clus02_RedType_GenTable_classic_FDR005.csv",header=T) %>% 
  select(GO.ID) %>% 
  slice(1:3) %>% 
  mutate("Clus"="RC02")
Clus03_Red <- read.csv("kmeans/ClusterInfo/topGO/Clus03_RedType_GenTable_classic_FDR005.csv",header=T) %>% 
  select(GO.ID) %>% 
  slice(1:3) %>% 
  mutate("Clus"="RC03")
Clus04_Red <- read.csv("kmeans/ClusterInfo/topGO/Clus04_RedType_GenTable_classic_FDR005.csv",header=T) %>% 
  select(GO.ID) %>% 
  slice(1:3) %>% 
  mutate("Clus"="RC04")
Clus05_Red <- read.csv("kmeans/ClusterInfo/topGO/Clus05_RedType_GenTable_classic_FDR005.csv",header=T) %>% 
  select(GO.ID) %>% 
  slice(1:3) %>% 
  mutate("Clus"="RC05")
Clus06_Red <- read.csv("kmeans/ClusterInfo/topGO/Clus06_RedType_GenTable_classic_FDR005.csv",header=T) %>% 
  select(GO.ID) %>% 
  slice(1:3) %>% 
  mutate("Clus"="RC06")

#Green
Clus01_Gre <- read.csv("kmeans/ClusterInfo/topGO/Clus01_GreenType_GenTable_classic_FDR005.csv",header=T) %>% 
  select(GO.ID) %>% 
  slice(1:3) %>% 
  mutate("Clus"="GC01")
Clus02_Gre <- read.csv("kmeans/ClusterInfo/topGO/Clus02_GreenType_GenTable_classic_FDR005.csv",header=T) %>% 
  select(GO.ID) %>% 
  slice(1:3) %>% 
  mutate("Clus"="GC02")
Clus03_Gre <- read.csv("kmeans/ClusterInfo/topGO/Clus03_GreenType_GenTable_classic_FDR005.csv",header=T) %>% 
  select(GO.ID) %>% 
  slice(1:3) %>% 
  mutate("Clus"="GC03")
Clus04_Gre <- read.csv("kmeans/ClusterInfo/topGO/Clus04_GreenType_GenTable_classic_FDR005.csv",header=T) %>% 
  select(GO.ID) %>% 
  slice(1:3) %>% 
  mutate("Clus"="GC04")
Clus05_Gre <- read.csv("kmeans/ClusterInfo/topGO/Clus05_GreenType_GenTable_classic_FDR005.csv",header=T) %>% 
  select(GO.ID) %>% 
  slice(1:3) %>% 
  mutate("Clus"="GC05")
Clus06_Gre <- read.csv("kmeans/ClusterInfo/topGO/Clus06_GreenType_GenTable_classic_FDR005.csv",header=T) %>% 
  select(GO.ID) %>% 
  slice(1:3) %>% 
  mutate("Clus"="GC06")

#Combine_Cluster
#Ver. select cluster
GO<- rbind.data.frame(Clus01_Red, Clus02_Red, Clus03_Red, Clus04_Red, Clus05_Red, Clus06_Red, 
                      Clus01_Gre, Clus02_Gre, Clus03_Gre, Clus04_Gre, Clus05_Gre, Clus06_Gre) 

#join MetaData
GO_Clus01_Red <- left_join(GO, 
                           read.csv("kmeans/ClusterInfo/topGO/Clus01_RedType_GenTable_classic.csv",header=T), 
                           by="GO.ID", na_matches="never")%>% 
  dplyr::select(GO.ID, Term, Pvalue) 

GO_Clus02_Red <- left_join(GO, 
                           read.csv("kmeans/ClusterInfo/topGO/Clus02_RedType_GenTable_classic.csv",header=T), 
                           by="GO.ID", na_matches="never")%>% 
  dplyr::select(GO.ID, Term, Pvalue) 

GO_Clus03_Red <- left_join(GO, 
                           read.csv("kmeans/ClusterInfo/topGO/Clus03_RedType_GenTable_classic.csv",header=T), 
                           by="GO.ID", na_matches="never")%>% 
  dplyr::select(GO.ID, Term, Pvalue) 

GO_Clus04_Red <- left_join(GO, 
                           read.csv("kmeans/ClusterInfo/topGO/Clus04_RedType_GenTable_classic.csv",header=T), 
                           by="GO.ID", na_matches="never")%>% 
  dplyr::select(GO.ID, Term, Pvalue) 

GO_Clus05_Red <- left_join(GO, 
                           read.csv("kmeans/ClusterInfo/topGO/Clus05_RedType_GenTable_classic.csv",header=T), 
                           by="GO.ID", na_matches="never")%>% 
  dplyr::select(GO.ID, Term, Pvalue) 

GO_Clus06_Red <- left_join(GO, 
                           read.csv("kmeans/ClusterInfo/topGO/Clus06_RedType_GenTable_classic.csv",header=T), 
                           by="GO.ID", na_matches="never")%>% 
  dplyr::select(GO.ID, Term, Pvalue) 

GO_Clus01_Gre <- left_join(GO, 
                           read.csv("kmeans/ClusterInfo/topGO/Clus01_GreenType_GenTable_classic.csv",header=T), 
                           by="GO.ID", na_matches="never")%>% 
  dplyr::select(GO.ID, Term, Pvalue) 

GO_Clus02_Gre <- left_join(GO, 
                           read.csv("kmeans/ClusterInfo/topGO/Clus02_GreenType_GenTable_classic.csv",header=T), 
                           by="GO.ID", na_matches="never")%>% 
  dplyr::select(GO.ID, Term, Pvalue) 

GO_Clus03_Gre <- left_join(GO, 
                           read.csv("kmeans/ClusterInfo/topGO/Clus03_GreenType_GenTable_classic.csv",header=T), 
                           by="GO.ID", na_matches="never")%>% 
  dplyr::select(GO.ID, Term, Pvalue) 

GO_Clus04_Gre <- left_join(GO, 
                           read.csv("kmeans/ClusterInfo/topGO/Clus04_GreenType_GenTable_classic.csv",header=T), 
                           by="GO.ID", na_matches="never")%>% 
  dplyr::select(GO.ID, Term, Pvalue) 

GO_Clus05_Gre <- left_join(GO, 
                           read.csv("kmeans/ClusterInfo/topGO/Clus05_GreenType_GenTable_classic.csv",header=T), 
                           by="GO.ID", na_matches="never")%>% 
  dplyr::select(GO.ID, Term, Pvalue) 

GO_Clus06_Gre <- left_join(GO, 
                           read.csv("kmeans/ClusterInfo/topGO/Clus06_GreenType_GenTable_classic.csv",header=T), 
                           by="GO.ID", na_matches="never")%>% 
  dplyr::select(GO.ID, Term, Pvalue) 


#Wide for heatmap
GO.wide <- cbind.data.frame(GO_Clus01_Red, 
                            GO_Clus02_Red %>% dplyr::select(Pvalue),
                            GO_Clus03_Red %>% dplyr::select(Pvalue),
                            GO_Clus04_Red %>% dplyr::select(Pvalue),
                            GO_Clus05_Red %>% dplyr::select(Pvalue),
                            GO_Clus06_Red %>% dplyr::select(Pvalue),
                            GO_Clus01_Gre %>% dplyr::select(Pvalue),
                            GO_Clus02_Gre %>% dplyr::select(Pvalue),
                            GO_Clus03_Gre %>% dplyr::select(Pvalue),
                            GO_Clus04_Gre %>% dplyr::select(Pvalue),
                            GO_Clus05_Gre %>% dplyr::select(Pvalue),
                            GO_Clus06_Gre %>% dplyr::select(Pvalue)
) 


colnames(GO.wide) <- c("GO.ID","Term",
                       "RC01","RC02","RC03","RC04","RC05","RC06",
                       "GC01","GC02","GC03","GC04","GC05","GC06")
GO.wide <- GO.wide %>% dplyr::distinct(Term,.keep_all=TRUE) 

#NA = 0.1
GO.wide[is.na(GO.wide)] <- 0.1

#as.numeric
GO.wide <- GO.wide %>% as.data.frame()
for (i in 3:14){
  GO.wide[,i] <- as.numeric(GO.wide[,i])
}

#Matrix for heatmap
Mat <- as.matrix(-log10(GO.wide[3:ncol(GO.wide)]))
rownames(Mat) <- GO.wide$Term

#select col pallete
colpal <- colorRampPalette(c("white","yellow2","orangered"))#light blue,black,yellow

#Setting label for heatmap
#Select cluster
Light <- c(rep(1,6), rep(2,6))

#visualize by complexHeatmap
pdf(file =paste("kmeans/Fig/",print(today,quate=F),"_Heatmap_GOenrich_SelectCluster_topGO_classic_Top3.pdf", sep=""),width=3.85, height=4.0)

Heatmap(Mat,
        col=colpal(10),
        #col=brewr.pal(9, "Reds"),
        column_order=sort(colnames(Mat)),
        column_split=Light, 
        top_annotation = HeatmapAnnotation(foo=anno_block(gp=gpar(fill=c("deeppink","green")),
                                                          labels=c("Red-type","Green-type"),
                                                          labels_gp=gpar(col=c("white","black"), fontsize=8))),
        border = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=8),
        rect_gp = gpar(col="black", lwd=1),
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        show_heatmap_legend = FALSE
        #heatmap_legend_param = list(title="-log10(FDR)", title_position="lefttop-rot")
)

dev.off()

#Output legends
pdf(file =paste("kmeans/Fig/",print(today,quate=F),"_Heatmap_GOenrich_SelectCluster_topGO_classic_Top3_Label.pdf", sep=""),width=4, height=5.5)

Heatmap(Mat,
        col=colpal(10),
        #col=brewr.pal(9, "Reds"),
        column_order=sort(colnames(Mat)),
        column_split=Light, 
        top_annotation = HeatmapAnnotation(foo=anno_block(gp=gpar(fill=c("deeppink","green")),
                                                          labels=c("Red-type","Green-type"),
                                                          labels_gp=gpar(col=c("white","black"), fontsize=8))),
        border = TRUE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize=8),
        column_names_gp = gpar(fontsize=8),
        rect_gp = gpar(col="black", lwd=1),
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        #show_heatmap_legend = FALSE,
        heatmap_legend_param = list(title="-log10(FDR)", title_position="lefttop-rot")
)

dev.off()



