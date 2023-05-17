##Fig.6d

#Upset plot

#Current date
today <- Sys.Date()
print(today, quote=F)

MLtoday <- "XXXX-XX-XX"

#load packages
library(tidyverse)
library(UpSetR)
library(ggupset)
library(ComplexHeatmap)

#met data
pheno <- read.csv("Pheno/Light_MorPheno.csv", header=T, check.names = FALSE)

#make empty matrix
FinalVar.df<-matrix(NA,ncol=2,nrow=0) #or 11

#make the dataframe of final-selected gens
for (i in 4:ncol(pheno)){
  print(i)
  #select phenp
  pheno.sel <- pheno %>% select(1,3,i)
  
  Trait=names(pheno.sel)[3]
  
  FinalVar<-read.csv(paste(MLtoday,"_",Algo,"/",Trait,"/BorutaAnotation/70/FinalVar.csv", sep=""), row.names = 1) %>% 
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
FinalVar.df.wide.Ano <- left_join(FinalVar.df.wide, GenesAno, by="GeneID", na_matches="never") 
FinalVar.df.wide.Ano  <- FinalVar.df.wide.Ano %>% distinct(GeneID, .keep_all=TRUE)

#output
write.csv(FinalVar.df.wide.Ano, paste(MLtoday,"_",Algo,"/Table_UpsetPlot/Table_for_UpsetPlot_Pheno.csv", sep=""), row.names = FALSE)


Growth = FinalVar.df %>% filter(Trait=="Growth") %>% select(GeneID)
Leaf_Width = FinalVar.df %>% filter(Trait=="Leaf_Width") %>% select(GeneID)
Leaf_Length = FinalVar.df %>% filter(Trait=="Leaf_Length") %>% select(GeneID) 
RLW = FinalVar.df %>% filter(Trait=="RLW") %>% select(GeneID) 
Leaf_Number = FinalVar.df %>% filter(Trait=="Leaf_Number") %>% select(GeneID) 

lt <- list(Growth = Growth[,1],
           Leaf_Width = Leaf_Width[,1],
           Leaf_Length = Leaf_Length[,1],
           RLW = RLW[,1],
           Leaf_Number = Leaf_Number[,1]
           )


#make the combination matrix
m1 <-make_comb_mat(lt, mode="distinct") #distinct or intersect or union
m2 <-make_comb_mat(lt, mode="intersect")
set_name(m1)
set_size(m1)
comb_size(m1)
comb_degree(m1)

#make the plot
sets <- c("Growth", "Leaf_Width","Leaf_Length","RLW","Leaf_Number")

#distinct
#Fig.6d

pdf(paste(MLtoday,"_",Algo,"/Fig/",today,"_UpsetPlot_AllPheno_distinct.pdf", sep=""), width=5, height=4)

UpSet(m1,
      set_order = sets,
      comb_col=c("deeppink2","dodgerblue3")[comb_degree(m1)],
      top_annotation = upset_top_annotation(m1, add_numbers=TRUE, 
                                                 gp=gpar(col="black", 
                                                         fill=c("deeppink2","dodgerblue3")[comb_degree(m1)]),
                                            ),
      #right_annotation = upset_right_annotation(m1, add_numbers=TRUE),
      right_annotation = rowAnnotation(
        "Set size"=anno_barplot(set_size(m1),
                                border = FALSE,
                                gp = gpar(col="black", fill="grey35"),
                                width = unit(2, "cm")
                                )
      ))

dev.off()

