##Fig.5d

#Upset plot

#Current date
today <- Sys.Date()
print(today, quote=F)

#library
library(tidyverse)
library(ComplexHeatmap)

#data import
lt <- list(RC01 = read.csv("kmeans/ClusterInfo/Clus01_Exp_RedType.csv",header=TRUE)[,38],
           RC02 = read.csv("kmeans/ClusterInfo/Clus02_Exp_RedType.csv",header=TRUE)[,38], 
           RC03 = read.csv("kmeans/ClusterInfo/Clus03_Exp_RedType.csv",header=TRUE)[,38],
           RC04 = read.csv("kmeans/ClusterInfo/Clus04_Exp_RedType.csv",header=TRUE)[,38],
           RC05 = read.csv("kmeans/ClusterInfo/Clus05_Exp_RedType.csv",header=TRUE)[,38],
           RC06 = read.csv("kmeans/ClusterInfo/Clus06_Exp_RedType.csv",header=TRUE)[,38],
           GC01 = read.csv("kmeans/ClusterInfo/Clus01_Exp_GreenType.csv",header=TRUE)[,22], 
           GC02 = read.csv("kmeans/ClusterInfo/Clus02_Exp_GreenType.csv",header=TRUE)[,22], 
           GC03 = read.csv("kmeans/ClusterInfo/Clus03_Exp_GreenType.csv",header=TRUE)[,22], 
           GC04 = read.csv("kmeans/ClusterInfo/Clus04_Exp_GreenType.csv",header=TRUE)[,22], 
           GC05 = read.csv("kmeans/ClusterInfo/Clus05_Exp_GreenType.csv",header=TRUE)[,22], 
           GC06 = read.csv("kmeans/ClusterInfo/Clus06_Exp_GreenType.csv",header=TRUE)[,22] 
           )

#make the combination matrix
m1 <-make_comb_mat(lt, mode="distinct") #distinct or intersect or union
m2 <-make_comb_mat(lt, mode="intersect")
set_name(m1)
set_size(m1)
comb_size(m1)
comb_degree(m1)

#make the plot
sets <- c("Chicoric_acid", "Chlrogenic_acid","Cy3_3MG","Cy3_6MG","Cy3G","Q3_6MbGG","Q3_6MGG","Q3_6MG","Q3G")

#distinct
pdf(file="kmeans/Fig/Fig5d_UpsetPlot_kmeans_Cluster.pdf",width=5, height=4)

UpSet(m1,
      set_order=rownames(m1),
      comb_col="dodgerblue3",
      top_annotation = upset_top_annotation(m1, add_numbers=FALSE,gp=gpar(col="black", fill="dodgerblue3")),
      right_annotation = rowAnnotation(
        "Set size"=anno_barplot(set_size(m1),
                                border = FALSE,
                                gp = gpar(col="black", fill="grey35"),
                                width = unit(2, "cm")
                                )
      )
      )

dev.off()



