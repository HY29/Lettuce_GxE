##Fig.8b

#Line graph by ggplot2

#Load library
library(tidyverse)
library(plotly)
library(edgeR)
library(viridis)
library(RColorBrewer)

#Current data
today <- Sys.Date()
print(today, quate=F)

#Input expression data (TPM) matrix (RSEM_output)
count <- read.csv("GeneExp/Light_Cultivar/GeneExpressionMatrix.count.csv", header=T, row.names = 1)
#Gene annotation Inf
GenesAno<-read.csv("GeneInfo/Lsativa_467_v5.annotation_info.txt", header=T, sep="\t")
PhenylAno<-read.csv("GeneInfo/PhenylPropanoid_Pathway_Lettuce_Genes.csv", header=T)
FlavoAno<-read.csv("GeneInfo/Flavonoid_Pathway_Lettuce_Genes.csv", header=T)

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
TMMnorm.log <- bind_cols(rownames(TMMnorm.log) %>% as.data.frame(), TMMnorm.log)
colnames(TMMnorm.log)[1] <- "GeneID"

#join targeted gene
Phenyl.Exp <- left_join(PhenylAno, TMMnorm.log, by="GeneID", na_matches="never") %>% na.omit()
Flavo.Exp <- left_join(FlavoAno, TMMnorm.log, by="GeneID", na_matches="never") %>% na.omit()

#MwtaData
MetaData <- read.csv("MetaData/SampleLabel_RLL4.csv", header=T)

#long
Phenyl.Exp.long <- Phenyl.Exp %>% 
  pivot_longer(4:ncol(Phenyl.Exp), names_to="Sample",values_to="RelExp")

Flavo.Exp.long <- Flavo.Exp %>% 
  pivot_longer(4:ncol(Flavo.Exp), names_to="Sample",values_to="RelExp")

#join MetaData
Phenyl.Exp.long.Meta <- left_join(MetaData, Phenyl.Exp.long, by="Sample", na_matches="never") %>% na.omit()
Flavo.Exp.long.Meta <- left_join(MetaData, Flavo.Exp.long, by="Sample", na_matches="never") %>% na.omit()

#combine
Pathway.Exp <- rbind(Phenyl.Exp.long.Meta, Flavo.Exp.long.Meta) 
#"red" filter
Pathway.Exp <- Pathway.Exp %>% 
  filter(Color=="Red")
#Setting of group turn
Pathway.Exp$Light <- factor(Pathway.Exp$Light, levels=c("WhiteLED", "RedLED","BlueLED","FL"))
Pathway.Exp$GeneSymbol <- factor(Pathway.Exp$GeneSymbol, levels=c("PAL", "C4H","4CL","HCT","C3H","CAS","CHS", "CHI","F3H","F3'H","FLS","DFR","ANS","UFGT"))
Pathway.Exp$Genotype <- factor(Pathway.Exp$Genotype, levels=c("RLL4","rll4"))

#ggplot2
#Flavonoid
p <- ggplot(Pathway.Exp, aes(x=Light,y=RelExp, fill=Genotype))+
  guides(color=FALSE, fill=FALSE, shape=FALSE)+
  stat_boxplot(color="black", geom="errorbar", width=0.6, position = position_dodge(width=0.6))+
  geom_boxplot(color="black", width=0.5,position = position_dodge(width=0.6), outlier.colour = NA)+
  geom_jitter(color="black", shape=21, size=0.8, position = position_dodge(width=0.6))+
  facet_wrap(~GeneSymbol+GeneID, ncol=6, scales="free")+
  theme_bw()+
  theme(text=element_text(size=10),
        axis.text=element_text(size=10,color="black"),
        #axis.text.x = element_text(size=10, color="black", angle=90, hjust=1, vjust=0.5),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 10, face = "plain"),
        legend.text = element_text(size = 10, face = "italic"))+
  scale_color_brewer(palette="Set2")+
  scale_fill_brewer(palette="Set2")+
  labs(y="Relative expression")

#plot
print(p)

ggsave("Fig/Fig8b_BoxPlot_GeneExp_FlavonoidPathway_Light.pdf", width=11, height=8)


