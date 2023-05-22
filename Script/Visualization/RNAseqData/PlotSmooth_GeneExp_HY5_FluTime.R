##Fig.S11e

#Line graph by ggplot2

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
HY5<-read.csv("GeneInfo/HY5_Lettuce_Genes.csv", header=T)

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
HY5.Exp <- left_join(HY5, TMMnorm.log, by="GeneID", na_matches="never") %>% na.omit()

#MwtaData
MetaData <- read.csv("MetaData/SampleLabel_Time.csv", header=T)

#long
HY5.Exp.long <- HY5.Exp %>% 
  pivot_longer(3:ncol(HY5.Exp), names_to="Sample",values_to="RelExp")


#join MetaData
HY5.Exp.long.Meta <- left_join(MetaData, HY5.Exp.long, by="Sample", na_matches="never") %>% na.omit()

#Mean
HY5.Exp.long.Meta.Mean <- HY5.Exp.long.Meta %>% 
  group_by(Time, Cultivar, Color, GeneSymbol, GeneID) %>% 
  summarize(Mean=mean(RelExp), SD=sd(RelExp)) 

#Setting of group turn
HY5.Exp.long.Meta$Genotype <- factor(HY5.Exp.long.Meta$Genotype, levels=c("RLL4","rll4"))

#ggplot2
#HY5
p <- ggplot(HY5.Exp.long.Meta, aes(x=Time,y=RelExp, color=Cultivar, fill=Cultivar))+
  guides(color=FALSE, fill=FALSE)+
  geom_point(size=0.5, shape=1, alpha=0.5)+
  #geom_line(aes(group=Cultivar),size=1)+
  #stat_summary(fun=mean,geom="line",size=0.1,color="blue")+
  #geom_line(stat="identity",color="blue",size=0.5)+
  #geom_bar(stat="identity",fill="cyan4",color="black", width=0.7,position=position_dodge(0.6))+
  #geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD, width=5),position=position_dodge(0.6),size=0.5, width=0.4)+
  stat_smooth(method="loess", size=0.6)+
  facet_wrap(~GeneSymbol+GeneID, ncol=6, scales="free")+
  theme_bw()+
  theme(text=element_text(size=10),axis.text=element_text(size=10,color="black"))+
  scale_color_manual(values=c("violetred3","violet", "Green","green4"))+
  scale_fill_manual(values=c("violetred3","violet", "Green","green4"))+
  labs(y="Relative expression",x="Time after fluorescent light irradiation (h)")

#plot
print(p)

ggsave("Fig/FigS11e_PlotSmooth_GeneExp_HY5_FluTimeSeries.pdf", width=3.2, height=2.3)




