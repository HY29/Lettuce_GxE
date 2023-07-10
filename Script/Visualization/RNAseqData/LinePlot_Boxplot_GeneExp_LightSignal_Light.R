##Fig.S15a,b

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
count <- read.csv("GeneExp/Light_Cultivar/GeneExpressionMatrix.count.csv", header=T, row.names = 1, sep="\t")
#Gene annotation Inf
GenesAno<-read.csv("GeneInfo/Lsativa_467_v5.annotation_info.txt", header=T, sep="\t")
LightSig<-read.csv("GeneInfo/Light_signal_Lettuce_Genes.csv", header=T)

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
LightSig.Exp <- left_join(LightSig, TMMnorm.log, by="GeneID", na_matches="never") %>% na.omit()

#MwtaData
MetaData <- read.csv("MetaData/SampleLabel_RLL4.csv", header=T)

#long
LightSig.Exp.long <- LightSig.Exp %>% 
  pivot_longer(3:ncol(LightSig.Exp), names_to="Sample",values_to="RelExp")


#join MetaData
LightSig.Exp.long.Meta <- left_join(MetaData, LightSig.Exp.long, by="Sample", na_matches="never") %>% na.omit()

#Mean
LightSig.Exp.long.Meta.Mean <- LightSig.Exp.long.Meta %>% 
  group_by(Light, Cultivar, Color, GeneName, GeneID) %>% 
  summarize(Mean=mean(RelExp), SD=sd(RelExp)) 

#Setting of group turn
LightSig.Exp.long.Meta.Mean$Light <- factor(LightSig.Exp.long.Meta.Mean$Light, levels=c("WhiteLED", "RedLED","BlueLED","FL"))
LightSig.Exp.long.Meta$Light <- factor(LightSig.Exp.long.Meta$Light, levels=c("WhiteLED", "RedLED","BlueLED","FL"))
LightSig.Exp.long.Meta$GeneName <- factor(LightSig.Exp.long.Meta$GeneName, levels=c("CRY1","CRY2","UVR8","COP1", "SPA1","SPA2","RUP (RLL4)","BIC"))
LightSig.Exp.long.Meta$Genotype <- factor(LightSig.Exp.long.Meta$Genotype, levels=c("RLL4","rll4"))


#ggplot2
#Line plot

p1 <- ggplot(LightSig.Exp.long.Meta.Mean, aes(x=Light,y=Mean, color=Color, fill=Color, shape=Cultivar, group=Cultivar))+
  guides(color=FALSE, fill=FALSE, shape=FALSE)+
  geom_line(size=0.5, alpha=0.8)+
  geom_point(size=2, alpha=0.8)+
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), position=position_dodge(0),size=0.5, width=1, alpha=0.8)+
  facet_wrap(~GeneName+GeneID, ncol=4, scales="free")+
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

ggsave("Fig/FigS15a_LinePlot_GeneExp_LightSignal_Light.pdf", width=8, height=6)

#"red" filter
LightSig.Exp.long.Meta <- LightSig.Exp.long.Meta %>% 
  filter(Color=="Red")

##ggplot2
#boxplot
p2 <- ggplot(LightSig.Exp.long.Meta, aes(x=Light,y=RelExp, fill=Genotype))+
  guides(color=FALSE, fill=FALSE, shape=FALSE)+
  stat_boxplot(color="black", geom="errorbar", width=0.6, position = position_dodge(width=0.6))+
  geom_boxplot(color="black", width=0.5,position = position_dodge(width=0.6), outlier.colour = NA)+
  geom_jitter(color="black", shape=21, size=0.8, position = position_dodge(width=0.6))+
  facet_wrap(~GeneName+GeneID, ncol=4, scales="free")+
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

ggsave("Fig/FigS15b_BoxPlot_GeneExp_LightSig_Light.pdf", width=8, height=6)



