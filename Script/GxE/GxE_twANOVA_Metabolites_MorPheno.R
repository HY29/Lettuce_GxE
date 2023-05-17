#########################
## GxE analysis
## two-way ANOVA
#########################

##Fig.S4

#Current date
today <- Sys.Date()
print(today, quote=F)

#load packages
library(tidyverse)
library(ggrepel)

#data import
pheno <- read.csv("Pheno/Light_MorPheno.csv",header=TRUE, check.names=FALSE) %>% 
  na.omit()
met <- read.csv("Pheno/Light_Metabolites.csv",header=TRUE, check.names=FALSE)

#make matrix
p.mat.pheno <- matrix(NA, ncol=3, nrow=ncol(pheno)-3)
p.mat.met <- matrix(NA, ncol=3, nrow=ncol(met)-4)
colnames(p.mat.pheno) <- c("Genotypes","Light","Interaction")
colnames(p.mat.met) <- c("Genotypes","Light","Interaction")

#pheno
for (i in 4:ncol(pheno)){
  print(i)
  Trait=names(pheno)[i]
  
  twANOVA <- aov(pheno[,i]~Cultivar*Light, data=pheno)
  summary(twANOVA)
  
  p.mat.pheno[(i-3),1] <- summary(twANOVA)[[1]][["Pr(>F)"]][1]
  p.mat.pheno[(i-3),2] <- summary(twANOVA)[[1]][["Pr(>F)"]][2]
  p.mat.pheno[(i-3),3] <- summary(twANOVA)[[1]][["Pr(>F)"]][3]
}

#mat
for (i in 5:ncol(met)){
  print(i)
  Trait=names(met)[i]
  
  twANOVA <- aov(met[,i]~Cultivar*Light, data=met)
  summary(twANOVA)
  
  p.mat.met[(i-4),1] <- summary(twANOVA)[[1]][["Pr(>F)"]][1]
  p.mat.met[(i-4),2] <- summary(twANOVA)[[1]][["Pr(>F)"]][2]
  p.mat.met[(i-4),3] <- summary(twANOVA)[[1]][["Pr(>F)"]][3]
}

p.mat.pheno <- bind_cols(colnames(pheno[4:ncol(pheno)]), p.mat.pheno %>% as.data.frame())
colnames(p.mat.pheno)[1]<-"Trait"

p.mat.met <- bind_cols(colnames(met[5:ncol(met)]), p.mat.met %>% as.data.frame())
colnames(p.mat.met)[1]<-"Trait"

p.mat <- rbind.data.frame(p.mat.met, p.mat.pheno)

p.mat.long <- p.mat %>% 
  pivot_longer(2:4, names_to="Stat", values_to="Pvalue")

p.mat.adj <- bind_cols(p.mat.long, p.adjust(p.mat.long$Pvalue, method="fdr"))
colnames(p.mat.adj)[4]<-"FDR"
p.mat.adj <- bind_cols(p.mat.adj, -log10(p.mat.adj$FDR))
colnames(p.mat.adj)[5]<-"Minus_log10_FDR"

p.mat.adj$Stat<- factor(p.mat.adj$Stat, levels=c("Genotypes","Light","Interaction"))
                                         
#Visualization
#ggplot2 
p <- ggplot(p.mat.adj, aes(y=Minus_log10_FDR, x=Stat, fill=Stat))+
  #guides(color=FALSE, fill=FALSE, shape=FALSE)+
  geom_bar(stat = "identity")+
  geom_hline(aes(yintercept=-log10(0.05)),color="blue", linetype="solid",size=0.5)+ #solid, dashed, dotted
  #facet_grid(Cultivar~metabolites, scales = "free")+
  facet_wrap(~ Trait, ncol=5, scales="free")+
  theme_bw()+
  theme(text = element_text(size = 12), axis.text = element_text(size=12, color="black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=12, color="black", angle=90, hjust=1, vjust=0.5))+
  scale_color_brewer(palette="Set2")+
  scale_fill_brewer(palette="Set2")+
  #scale_color_manual(values=c("plum1","indianred4"))+
  #scale_fill_manual(values=c("plum1","indianred4"))+
  #scale_fill_viridis_d(option="cividis")+
  labs(y="-log10(FDR)")

#plot
print(p)

ggsave("Fig/FigS4_TwoWayANOVA_Light_Metabolites.pdf", width=8.5, height=6.5)


