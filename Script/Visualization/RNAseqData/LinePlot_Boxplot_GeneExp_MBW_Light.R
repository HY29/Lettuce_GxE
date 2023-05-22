##Fig.S12a,b

#Line graph & Boxplot by ggplot2

#Load library
library(tidyverse)
library(edgeR)
library(viridis)
library(RColorBrewer)

#Current data
today <- Sys.Date()
print(today, quate=F)

#Input expression data (TPM) matrix (RSEM_output)
count <- read.csv("GeneExp/Light_Cultivar/GeneExpressionMatrix.count.csv", header=T, row.names = 1, sep="")
#Gene annotation Inf
GenesAno<-read.csv("GeneInfo/Lsativa_467_v5.annotation_info.txt", header=T, sep="\t")
MBW<-read.csv("GeneInfo/MBW_complex_v2.csv", header=T)

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
MBW.Exp <- left_join(MBW, TMMnorm.log, by="GeneID", na_matches="never") %>% na.omit()

#MwtaData
MetaData <- read.csv("MetaData/SampleLabel_RLL4.csv", header=T)

#long
MBW.Exp.long <- MBW.Exp %>% 
  pivot_longer(6:ncol(MBW.Exp), names_to="Sample",values_to="RelExp")


#join MetaData
MBW.Exp.long.Meta <- left_join(MetaData, MBW.Exp.long, by="Sample", na_matches="never") %>% na.omit()

#Setting of group turn
MBW.Exp.long.Meta.Mean$Light <- factor(MBW.Exp.long.Meta.Mean$Light, levels=c("WhiteLED", "RedLED","BlueLED","FL"))
MBW.Exp.long.Meta$Light <- factor(MBW.Exp.long.Meta$Light, levels=c("WhiteLED", "RedLED","BlueLED","FL"))
MBW.Exp.long.Meta <- mutate(MBW.Exp.long.Meta, GeneSymbol=fct_inorder(GeneSymbol))
MBW.Exp.long.Meta <- mutate(MBW.Exp.long.Meta, GeneID=fct_inorder(GeneID))

#HY5.Exp.long.Meta$GeneSymbol <- factor(HY5.Exp.long.Meta$GeneSymbol, levels=c("PAL", "C4H","4CL","HCT","C3H","CAS","CHS", "CHI","F3H","F3'H","FLS","DFR","ANS","UFGT"))
MBW.Exp.long.Meta$Genotype <- factor(MBW.Exp.long.Meta$Genotype, levels=c("RLL4","rll4"))

#Mean
MBW.Exp.long.Meta.Mean <- MBW.Exp.long.Meta %>% 
  group_by(Light, Cultivar, Color, GeneSymbol, GeneID) %>% 
  summarize(Mean=mean(RelExp), SD=sd(RelExp)) 

#ggplot2
#Line plot
p1 <- ggplot(MBW.Exp.long.Meta.Mean, aes(x=Light,y=Mean, color=Color, fill=Color, shape=Cultivar, group=Cultivar))+
  guides(color=FALSE, fill=FALSE, shape=FALSE)+
  geom_line(size=0.5, alpha=0.8)+
  geom_point(size=2, alpha=0.8)+
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), position=position_dodge(0),size=0.5, width=1, alpha=0.8)+
  facet_wrap(~GeneSymbol+GeneID, ncol=6, scales="free")+
  theme_bw()+
  theme(text=element_text(size=10),
        axis.text=element_text(size=10,color="black"),
        #axis.text.x = element_text(size=10, color="black", angle=90, hjust=1, vjust=0.5),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())+
  scale_shape_manual(values=c(0,1,2,3,4,5,6,7,8,9,15,16,17,18))+
  scale_color_manual(values=c("green","deeppink2"))+
  labs(y="Relative expression")

#ggplotly
ggplotly(p1)

#plot
print(p1)

ggsave("Fig/FigS12a_LinePlot_GeneExp_MBW_Light.pdf", width=10, height=3.5)

#"red" filter
MBW.Exp.long.Meta <- MBW.Exp.long.Meta %>% 
  filter(Color=="Red")

#Mean
MBW.Exp.long.Meta.Mean <- MBW.Exp.long.Meta %>% 
  group_by(Light, Cultivar, Color, GeneSymbol, GeneID) %>% 
  summarize(Mean=mean(RelExp), SD=sd(RelExp)) 


##ggplot2
#boxplot
p2 <- ggplot(MBW.Exp.long.Meta, aes(x=Light,y=RelExp, fill=Genotype))+
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
print(p2)

ggsave("Fig/FigS12b_BoxPlot_GeneExp_MBW_Light.pdf", width=10, height=3.5)




