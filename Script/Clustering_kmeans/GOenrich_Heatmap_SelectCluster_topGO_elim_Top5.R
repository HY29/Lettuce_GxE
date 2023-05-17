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
Clus01_Red <- read.csv("kmeans/topGO/Clus01_RedType_GenTable_elim_FDR005.csv",header=T) %>% 
  dplyr::select(GO.ID,Term, Pvalue) %>% 
  dplyr::slice(1:5) %>% 
  dplyr::mutate("Clus"="RC01")
Clus02_Red <- read.csv("kmeans/topGO/Clus02_RedType_GenTable_elim_FDR005.csv",header=T) %>% 
  dplyr::select(GO.ID,Term, Pvalue) %>% 
  dplyr::slice(1:5) %>% 
  dplyr::mutate("Clus"="RC02")
Clus03_Red <- read.csv("kmeans/topGO/Clus03_RedType_GenTable_elim_FDR005.csv",header=T) %>% 
  dplyr::select(GO.ID,Term, Pvalue) %>% 
  dplyr::slice(1:5) %>% 
  dplyr::mutate("Clus"="RC03")
Clus04_Red <- read.csv("kmeans/topGO/Clus04_RedType_GenTable_elim_FDR005.csv",header=T) %>% 
  dplyr::select(GO.ID,Term, Pvalue) %>% 
  dplyr::slice(1:5) %>% 
  dplyr::mutate("Clus"="RC04")
Clus05_Red <- read.csv("kmeans/topGO/Clus05_RedType_GenTable_elim_FDR005.csv",header=T) %>% 
  dplyr::select(GO.ID,Term, Pvalue) %>% 
  dplyr::slice(1:5) %>% 
  dplyr::mutate("Clus"="RC05")
Clus06_Red <- read.csv("kmeans/topGO/Clus06_RedType_GenTable_elim_FDR005.csv",header=T) %>% 
  dplyr::select(GO.ID,Term, Pvalue) %>% 
  dplyr::slice(1:5) %>% 
  dplyr::mutate("Clus"="RC06")

#Green
Clus01_Gre <- read.csv("kmeans/topGO/Clus01_GreenType_GenTable_elim_FDR005.csv",header=T) %>% 
  dplyr::select(GO.ID,Term, Pvalue) %>% 
  dplyr::slice(1:5) %>% 
  dplyr::mutate("Clus"="GC01")
Clus02_Gre <- read.csv("kmeans/topGO/Clus02_GreenType_GenTable_elim_FDR005.csv",header=T) %>% 
  dplyr::select(GO.ID,Term, Pvalue) %>% 
  dplyr::slice(1:5) %>% 
  dplyr::mutate("Clus"="GC02")
Clus03_Gre <- read.csv("kmeans/topGO/Clus03_GreenType_GenTable_elim_FDR005.csv",header=T) %>% 
  dplyr::select(GO.ID,Term, Pvalue) %>% 
  dplyr::slice(1:5) %>% 
  dplyr::mutate("Clus"="GC03")
Clus04_Gre <- read.csv("kmeans/topGO/Clus04_GreenType_GenTable_elim_FDR005.csv",header=T) %>% 
  dplyr::select(GO.ID,Term, Pvalue) %>% 
  dplyr::slice(1:5) %>% 
  dplyr::mutate("Clus"="GC04")
Clus05_Gre <- read.csv("kmeans/topGO/Clus05_GreenType_GenTable_elim_FDR005.csv",header=T) %>% 
  dplyr::select(GO.ID,Term, Pvalue) %>% 
  dplyr::slice(1:5) %>% 
  dplyr::mutate("Clus"="GC05")
Clus06_Gre <- read.csv("kmeans/topGO/Clus06_GreenType_GenTable_elim_FDR005.csv",header=T) %>% 
  dplyr::select(GO.ID,Term, Pvalue) %>% 
  dplyr::slice(1:5) %>% 
  dplyr::mutate("Clus"="GC06")

#Combine_Cluster
#Ver. select cluster
GO_Clus <- bind_rows(Clus02_Red, Clus03_Red, Clus04_Red, 
                    Clus02_Gre, Clus03_Gre) 

#Wide for heatmap
GO_Clus.wide <- GO_Clus %>% 
  pivot_wider(names_from=Clus, values_from=Pvalue) %>% 
  mutate("GC06"=NA)

#NA = 0.1
GO_Clus.wide[is.na(GO_Clus.wide)] <- 0.1

#Matrix for heatmap
Mat <- as.matrix(-log10(GO_Clus.wide[3:ncol(GO_Clus.wide)]))
rownames(Mat) <- GO_Clus.wide$Term

#select col pallete
colpal <- colorRampPalette(c("white","yellow2","orangered"))#light blue,black,yellow

#Setting label for heatmap
#Select cluster
Light <- c(rep(1,3), rep(2,3))

#visualize by complexHeatmap
pdf(file = "kmeans/Fig/Fig5e_Heatmap_GOenrich_SelectCluster_topGO_elim_Top5.pdf",width=4.5, height=4)

Heatmap(Mat,
        col=colpal(10),
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
        #show_heatmap_legend = FALSE,
        heatmap_legend_param = list(title="-log10(FDR)", title_position="lefttop-rot")
        )
        

dev.off()

