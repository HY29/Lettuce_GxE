##Fig.S9a

#Line graph by ggplot2

#library
library(tidyverse)
library(edgeR)
library(viridis)
library(RColorBrewer)

#Current data
today <- Sys.Date()
print(today, quate=F)

#Input expression data (TPM) matrix (RSEM_output)
count <- read.csv("GeneExp/UVB/GeneExpressionMatrix.count.UVB.csv", header=T, row.names = 1)
#Gene annotation Inf
GenesAno<-read.csv("GeneInfo/Lsativa_467_v5.annotation_info.txt", header=T, sep="\t")
PhenylAno<-read.csv("GeneInfo/PhenylPropanoid_Pathway_Lettuce_Genes.csv", header=T)
FlavoAno<-read.csv("GeneInfo/Flavonoid_Pathway_Lettuce_Genes.csv", header=T)

#except Ls21v2-084 (low-depth read)
count <- count %>% dplyr::select(-Ls21v2.084.genes.results)

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
MetaData <- read.csv("MetaData/SampleLabel_UVB.csv", header=T)

#long
Phenyl.Exp.long <- Phenyl.Exp %>% 
  pivot_longer(4:ncol(Phenyl.Exp), names_to="Sample",values_to="RelExp")

Flavo.Exp.long <- Flavo.Exp %>% 
  pivot_longer(4:ncol(Flavo.Exp), names_to="Sample",values_to="RelExp")

#join MetaData
Phenyl.Exp.long.Meta <- left_join(MetaData, Phenyl.Exp.long, by="Sample", na_matches="never") %>% na.omit()
Flavo.Exp.long.Meta <- left_join(MetaData, Flavo.Exp.long, by="Sample", na_matches="never") %>% na.omit()

#Mean
Phenyl.Exp.long.Meta.Mean <- Phenyl.Exp.long.Meta %>% 
  group_by(Light, Cultivar, Color, GeneSymbol, GeneID) %>% 
  summarize(Mean=mean(RelExp), SD=sd(RelExp)) 

Flavo.Exp.long.Meta.Mean <- Flavo.Exp.long.Meta %>% 
  group_by(Light, Cultivar, Color, GeneSymbol, GeneID) %>% 
  summarize(Mean=mean(RelExp), SD=sd(RelExp)) 

#Setting of group turn
Phenyl.Exp.long.Meta.Mean$GeneSymbol <- factor(Phenyl.Exp.long.Meta.Mean$GeneSymbol, levels=c("PAL", "C4H","4CL","HCT","C3H","CAS"))
#Setting of group turn
Flavo.Exp.long.Meta.Mean$GeneSymbol <- factor(Flavo.Exp.long.Meta.Mean$GeneSymbol, levels=c("CHS", "CHI","F3H","F3'H","FLS","DFR","ANS","UFGT"))
#combine
Pathway.Exp <- rbind(Phenyl.Exp.long.Meta.Mean, Flavo.Exp.long.Meta.Mean)

#ggplot2
#Flavonoid pathway
p <- ggplot(Pathway.Exp, aes(x=Light,y=Mean, color=Color, fill=Color, shape=Cultivar, group=Cultivar))+
  geom_line(size=0.5, alpha=0.8)+
  geom_point(size=2, alpha=0.8)+
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), position=position_dodge(0),size=0.5, width=1, alpha=0.8)+
  facet_wrap(~GeneSymbol+GeneID, ncol=6, scales="free")+
  theme_bw()+
  theme(text=element_text(size=10),
        #axis.text=element_text(size=10,color="black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  scale_shape_manual(values=c(0:2,4:7,15))+
  scale_color_manual(values=c("green","deeppink2"))+
  #scale_fill_manual(values=c("Red","deeppink", "Green","Blue"))+
  labs(y="Relative expression")

#ggplotly
ggplotly(p)

#plot
print(p)

ggsave("Fig/FigS9a_LinePlot_GeneExp_FlavonoidPathway_UVB.pdf", width=11, height=8)
