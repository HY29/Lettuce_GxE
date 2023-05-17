##Fig.6c

#Upset plot

#Current date
today <- Sys.Date()
print(today,quote=F)

MLtoday <- "XXXX-XX-XX"

#load packages
library(tidyverse)
library(UpSetR)
library(ggupset)
library(ComplexHeatmap)

#met data
met <- read.csv("Pheno/Light_Metabolites.csv", header=T, check.names = FALSE)

#make empty matrix
FinalVar.df<-matrix(NA,ncol=2,nrow=0) #or 11

#make the dataframe of final-selected gens
for (i in 5:ncol(met.ori)){
  print(i)
  #select phenp
  pheno.sel <- met %>% select(1,4,i)
  
  Trait=names(pheno.sel)[3]
  
  FinalVar<-read.csv(paste(MLtoday,"_",Algo,"/",Trait,"/BorutaAnotation/100/FinalVar.csv", sep=""), row.names = 1) %>% 
    mutate(Trait=Trait)
  
  FinalVar.df <- rbind.data.frame(FinalVar.df, FinalVar)
} 

#for output
#Gene annotation Inf
GenesAno<-read.csv("GeneInfo/Lsativa_467_v5.annotation_info.txt", header=T, sep="\t")
PhenylAno<-read.csv("GeneInfo/PhenylPropanoid_Pathway_Lettuce_Genes.csv", header=T)
FlavoAno<-read.csv("GeneInfo/Flavonoid_Pathway_Lettuce_Genes.csv", header=T)
PathwayAno <- rbind(PhenylAno, FlavoAno)

#wide
FinalVar.df.wide <- FinalVar.df %>% pivot_wider(names_from = Trait, values_from=Trait)
#join targeted gene
FinalVar.df.wide.Ano <- left_join(FinalVar.df.wide, PathwayAno, by="GeneID", na_matches="never") 
#output
write.csv(FinalVar.df.wide.Ano, paste(MLtoday,"_",Algo,"/Table_UpsetPlot/Table_for_UpsetPlot_Met.csv", sep=""), row.names = FALSE)


Chicoric_acid = FinalVar.df %>% filter(Trait=="Chicoric acid") %>% select(GeneID)
Chlrogenic_acid = FinalVar.df %>% filter(Trait=="Chlrogenic acid") %>% select(GeneID)
Cy_3MG = FinalVar.df %>% filter(Trait=="Cy3_3MG") %>% select(GeneID) 
Cy_6MG = FinalVar.df %>% filter(Trait=="Cy3_6MG") %>% select(GeneID) 
Cy3G = FinalVar.df %>% filter(Trait=="Cy3G") %>% select(GeneID) 
Q3_6MbGG = FinalVar.df %>% filter(Trait=="Q3_6MbGG") %>% select(GeneID) 
Q3_6MGG = FinalVar.df %>% filter(Trait=="Q3_6MGG") %>% select(GeneID) 
Q3_6MG = FinalVar.df %>% filter(Trait=="Q3_6MG") %>% select(GeneID) 
Q3G = FinalVar.df %>% filter(Trait=="Q3G") %>% select(GeneID) 

lt <- list(Chicoric_acid = Chicoric_acid[,1],
           Chlrogenic_acid = Chlrogenic_acid[,1],
           Cy3_3MG = Cy_3MG[,1],
           Cy3_6MG = Cy_6MG[,1],
           Cy3G = Cy3G[,1],
           Q3_6MbGG = Q3_6MbGG[,1],
           Q3_6MGG = Q3_6MGG[,1],
           Q3_6MG = Q3_6MG[,1],
           Q3G = Q3G[,1])


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
#Fig.6c
pdf(paste(MLtoday,"_",Algo,"/Fig/",today,"_UpsetPlot_AllMet_distinct.pdf", sep=""), width=7, height=4)

UpSet(m1,
      set_order = sets,
      comb_col=c("deeppink2","dodgerblue3","tomato1","orange","springgreen2")[comb_degree(m1)],
      top_annotation = upset_top_annotation(m1, add_numbers=TRUE, 
                                                 gp=gpar(col="black", 
                                                         fill=c("deeppink2","dodgerblue3","tomato1","orange","springgreen2")[comb_degree(m1)]),
                                            ),
      #right_annotation = upset_right_annotation(m1, add_numbers=TRUE),
      right_annotation = rowAnnotation(
        "Set size"=anno_barplot(set_size(m1),
                                border = FALSE,
                                gp = gpar(col="black", fill="grey35"),
                                width = unit(2, "cm")
             ),
             Group=c("Chicoric_acid","Chlrogenic_acid",
                     "Cyanidin","Cyanidin","Cyanidin",
                     "Quercetin","Quercetin","Quercetin","Quercetin"))
      )

dev.off()

