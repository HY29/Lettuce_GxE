##Fig.S5a

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
Met <- read.csv("Pheno/FluTime_Metabolites.csv", header=T, check.names = FALSE)

#long
Met.long <- Met %>% 
  pivot_longer(7:ncol(Met), names_to="Metabolite",values_to="Value")

#ggplot2
ggplot(Met.long, aes(x=Time,y=Value, color=Cultivar, fill=Cultivar))+
  #guides(color=FALSE, fill=FALSE)+
  geom_point(size=0.5, shape=1, alpha=0.5)+
  stat_smooth(method="loess", size=0.6)+
  facet_wrap(~Metabolite, ncol=3, scales="free")+
  theme_bw()+
  theme(text=element_text(size=10),axis.text=element_text(size=10,color="black"))+
  scale_color_manual(values=c("violetred3","violet", "Green","green4"))+
  scale_fill_manual(values=c("violetred3","violet", "Green","green4"))+
  labs(y="mg/100 gFW",x="Time after fluorescent light irradiation (h)")

ggsave("Fig/FigS5a_PlotSmooth_Metabolites_FluTimeSeries.pdf", width=7, height=5)



