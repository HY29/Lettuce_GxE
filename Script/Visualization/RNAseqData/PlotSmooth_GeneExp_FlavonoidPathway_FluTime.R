##Fig.S5b

#Plot Smooth by ggplot2

#Load library
library(tidyverse)
library(edgeR)
library(viridis)
library(RColorBrewer)

#Current data
today <- Sys.Date()
print(today, quate=F)

#Input expression data (TPM) matrix (RSEM_output)
count <- read.csv("GeneExp/Flu_TimeSeries/GeneExpressionMatrix.count.Time.csv", header=T, row.names = 1)
#Gene annotation Inf
GenesAno<-read.csv("GeneInfo/Lsativa_467_v5.annotation_info.txt", header=T, sep="\t")
PhenylAno<-read.csv("GeneInfo/PhenylPropanoid_Pathway_Lettuce_Genes.csv", header=T)
FlavoAno<-read.csv("GeneInfo/Flavonoid_Pathway_Lettuce_Genes.csv", header=T)

#except Ls21v2-084 (low-depth read)
#count <- count %>% select(-Ls21v2.084.genes.results)

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
MetaData <- read.csv("MetaData/SampleLabel_Time.csv", header=T)

#long
Phenyl.Exp.long <- Phenyl.Exp %>% 
  pivot_longer(4:ncol(Phenyl.Exp), names_to="Sample",values_to="RelExp")

Flavo.Exp.long <- Flavo.Exp %>% 
  pivot_longer(4:ncol(Flavo.Exp), names_to="Sample",values_to="RelExp")

#join MetaData
Phenyl.Exp.long.Meta <- left_join(MetaData, Phenyl.Exp.long, by="Sample", na_matches="never") %>% na.omit()
Flavo.Exp.long.Meta <- left_join(MetaData, Flavo.Exp.long, by="Sample", na_matches="never") %>% na.omit()
#combine
Pathway.Exp <- rbind(Phenyl.Exp.long.Meta, Flavo.Exp.long.Meta )

#Setting of group turn
Pathway.Exp$GeneSymbol <- factor(Pathway.Exp$GeneSymbol, levels=c("PAL", "C4H","4CL","HCT","C3H","CAS","CHS", "CHI","F3H","F3'H","FLS","DFR","ANS","UFGT"))

#ggplot2
#PhenylPropanoid
ggplot(Pathway.Exp, aes(x=Time,y=RelExp, color=Cultivar, fill=Cultivar))+
  #guides(color=FALSE, fill=FALSE)+
  geom_point(size=0.5, shape=1, alpha=0.5)+
  stat_smooth(method="loess", size=0.6)+
  facet_wrap(~GeneSymbol+GeneID, ncol=6, scales="free")+
  theme_bw()+
  theme(text=element_text(size=10),axis.text=element_text(size=10,color="black"))+
  scale_color_manual(values=c("violetred3","violet", "Green","green4"))+
  scale_fill_manual(values=c("violetred3","violet", "Green","green4"))+
  labs(y="Relative expression",x="Time after fluorescent light irradiation (h)")
ggsave("Fig/FigS5b_PlotSmooth_GeneExp_FlavonoidPathway_FluTimeSeries.pdf", width=10, height=7)



