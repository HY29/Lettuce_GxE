##Fig.S13e

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
MetaData <- read.csv("MetaData/SampleLabel_Time.csv", header=T)

#long
LightSig.Exp.long <- LightSig.Exp %>% 
  pivot_longer(6:ncol(LightSig.Exp), names_to="Sample",values_to="RelExp")


#join MetaData
LightSig.Exp.long.Meta <- left_join(MetaData, LightSig.Exp.long, by="Sample", na_matches="never") %>% na.omit()
LightSig.Exp.long.Meta$GeneName <- factor(LightSig.Exp.long.Meta$GeneName, levels=c("CRY1","CRY2","UVR8","COP1", "SPA1","SPA2","RUP (RLL4)","BIC"))


#ggplot2
p1 <- ggplot(LightSig.Exp.long.Meta, aes(x=Time,y=RelExp, color=Cultivar, fill=Cultivar))+
  guides(color=FALSE, fill=FALSE)+
  geom_point(size=0.5, shape=1, alpha=0.5)+
  stat_smooth(method="loess", size=0.6)+
  facet_wrap(~GeneName+GeneID, ncol=6, scales="free")+
  theme_bw()+
  theme(text=element_text(size=10),axis.text=element_text(size=10,color="black"))+
  scale_color_manual(values=c("violetred3","violet", "Green","green4"))+
  scale_fill_manual(values=c("violetred3","violet", "Green","green4"))+
  labs(y="Relative expression",x="Time after fluorescent light irradiation (h)")

#plot 
print(p1)

ggsave("Fig/FigS13e_PlotSmooth_GeneExp_LightSig_FluTimeSeries.pdf", width=10.5, height=5)




