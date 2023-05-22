##Fig.S12e

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
MetaData <- read.csv("MetaData/SampleLabel_Time.csv", header=T)

#long
MBW.Exp.long <- MBW.Exp %>% 
  pivot_longer(6:ncol(MBW.Exp), names_to="Sample",values_to="RelExp")


#join MetaData
MBW.Exp.long.Meta <- left_join(MetaData, MBW.Exp.long, by="Sample", na_matches="never") %>% na.omit()

#Setting of group turn
MBW.Exp.long.Meta <- mutate(MBW.Exp.long.Meta, GeneSymbol=fct_inorder(GeneSymbol))
MBW.Exp.long.Meta <- mutate(MBW.Exp.long.Meta, GeneID=fct_inorder(GeneID))

#ggplot2
#HY5
p <- ggplot(MBW.Exp.long.Meta, aes(x=Time,y=RelExp, color=Cultivar, fill=Cultivar))+
  guides(color=FALSE, fill=FALSE)+
  geom_point(size=0.5, shape=1, alpha=0.5)+
  stat_smooth(method="loess", size=0.6)+
  facet_wrap(~GeneSymbol+GeneID, ncol=5, scales="free")+
  theme_bw()+
  theme(text=element_text(size=10),axis.text=element_text(size=10,color="black"))+
  scale_color_manual(values=c("violetred3","violet", "Green","green4"))+
  scale_fill_manual(values=c("violetred3","violet", "Green","green4"))+
  labs(y="Relative expression",x="Time after fluorescent light irradiation (h)")

#plot
print(p)

ggsave("Fig/FigS12e_PlotSmooth_GeneExp_MBW_FluTimeSeries.pdf", width=9, height=3.4)




