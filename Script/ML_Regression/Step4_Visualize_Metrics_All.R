#########################################################
### Calculation of R and RMSE value in genomic prediction
### Draft script 200811 HY
### Update 210217 HY
#########################################################

##Fig.6a, Fig.6b

#library
library(Metrics)
library(tidyverse)

##SetParameters
#Set seed
set.seed(123)
#Current date
today <- Sys.Date()
MLtoday <- "XXXX-XX-XX"
print(today,quote=F)
#Algorithms
Algo=print("RF_Boruta")
#Number of repetition
Repeat=100 

## Modelling based on selected feature----
#Pheno data
met <- read.csv("Pheno/Light_Metabolites.csv", header=T, check.names = FALSE)
pheno <- read.csv("Pheno/Light_MorPheno.csv", header=T, check.names = FALSE)

#empty matrix for R2.df
R2.df <- matrix(NA,nrow=0, ncol=6)
#make R2.df in all rep.
for (ThrVal in c(seq(10,100,10))){
  R2<-read.csv(paste(MLtoday,"_",Algo,"/PerformanceRawData/",ThrVal,"/AllMet_R2.csv", sep=""), header=T, check.names = FALSE) 
  colnames(R2)<-c("RepNo",colnames(met)[5:ncol(met)])
  R2.long <- R2 %>% 
  pivot_longer(2:ncol(R2), names_to = "Trait", values_to = "R2") %>% 
  mutate(ThrVal=ThrVal)
  
  R2.df <- rbind.data.frame(R2.df, R2.long)
}
#empty matrix for RMSE.df
RMSE.df <- matrix(NA,nrow=0, ncol=6)
#make RMSE.df in all rep.
for (ThrVal in c(seq(10,100,10))){
  RMSE<-read.csv(paste(MLtoday,"_",Algo,"/PerformanceRawData/",ThrVal,"/AllMet_RMSE.csv", sep=""), header=T, check.names = FALSE) 
  colnames(RMSE)<-c("RepNo",colnames(met)[5:ncol(met)])
  RMSE.long <- RMSE %>% 
    pivot_longer(2:ncol(RMSE), names_to = "Trait", values_to = "RMSE") %>% 
    mutate(ThrVal=ThrVal)
  
  RMSE.df <- rbind.data.frame(RMSE.df, RMSE.long)
}
#empty matrix for MAE.df
MAE.df <- matrix(NA,nrow=0, ncol=6)
#make MAE.df in all rep.
for (ThrVal in c(seq(10,100,10))){
  MAE<-read.csv(paste(MLtoday,"_",Algo,"/PerformanceRawData/",ThrVal,"/AllMet_MAE.csv", sep=""), header=T, check.names = FALSE) 
  colnames(MAE)<-c("RepNo",colnames(met)[5:ncol(met)])
  MAE.long <- MAE %>% 
    pivot_longer(2:ncol(MAE), names_to = "Trait", values_to = "MAE") %>% 
    mutate(ThrVal=ThrVal)
  
  MAE.df <- rbind.data.frame(MAE.df, MAE.long)
}

#join
Perform.met.df <- left_join(R2.df, RMSE.df, by=c("RepNo","Trait","ThrVal"), keep=FALSE) 
Perform.met.df <- left_join(Perform.met.df, MAE.df, by=c("RepNo","Trait","ThrVal"), keep=FALSE) 
#change col
Perform.met.df <- Perform.met.df[,c(1,5,3,6,2,4)]

#Perform pheno data
Perform.pheno.df <- read.csv("2022-05-25_RF_Boruta/PerformanceRawData/Pheno/AllMetrics_Pheno.csv", header=T, row.names=1,check.names = FALSE)
#Combine met & pheno df
Combine.df <- rbind(Perform.met.df, Perform.pheno.df)

#Setting of group turn
Combine.df$Trait <- factor(Combine.df$Trait, levels=c("Chicoric acid", "Chlrogenic acid","Cy3_3MG",
                                                      "Cy3_6MG","Cy3G","Q3_6MbGG",
                                                      "Q3_6MGG","Q3_6MG","Q3G",
                                                      "Growth", "Leaf_Length","Leaf_Number","Leaf_Width","RLW"))

#visualization of R2
#ggplot
#Fig.6a
ggplot(Combine.df, aes(x=ThrVal, y=R2))+
  geom_point(size=1, shape=1, color="dodgerblue2")+
  stat_smooth(method="loess",color="orangered1",size=0.6)+
  facet_wrap(~Trait, ncol=7)+
  theme_bw()+
  theme(text=element_text(size=12), 
        axis.text=element_text(size=12, color="black"),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        #axis.title.x = element_blank(),
        )+
  scale_x_continuous(breaks=seq(10,100,10), limits=c(0,110), expand=c(0,0))+
  labs(y=expression(paste("Prediction accuracy"," (",{R[2]},")", sep="")), x="Number of times selected in 100 repetitions")
ggsave(paste(MLtoday,"_",Algo,"/Fig/",today,"_Plot_R2_All.pdf",sep=""), width=10, height=4)

#Extract a number of FinalVal genes
#met data
met.ori <- read.csv("Pheno/21_lettuce_02_metabolites.csv", header=T, check.names = FALSE)

num.gene.mat<-matrix(NA,ncol=ncol(met.ori[,4:12]),nrow=1)
num.gene.mat.df<-matrix(NA,ncol=ncol(met.ori[,4:12]),nrow=0) #or 11

#rownames(sel.genenum.mat) <- c(1, seq(10,100,10))
#rownames(num.gene.mat) <- c(seq(10,100,10))
colnames(num.gene.mat) <- colnames(met.ori[,4:12])
colnames(num.gene.mat.df) <- colnames(met.ori[,4:12])

#Extract each gene number
for (ThrVal in c(seq(10,100,10))){
  print(ThrVal)
for (i in 4:ncol(met.ori)){
  print(paste(ThrVal,i))
  #select phenp
  pheno.sel <- met.ori %>% select(1,3,i)
  
  Trait=names(pheno.sel)[3]
  
  FinalVar<-read.csv(paste(MLtoday,"_",Algo,"/",Trait,"/BorutaAnotation/",ThrVal,"/FinalVar.csv", sep=""), row.names = 1)
  num.gene.mat[1,i-3] <- nrow(FinalVar)
  } 
  num.gene.mat.df <- rbind.data.frame(num.gene.mat.df, num.gene.mat)
}
colnames(num.gene.mat.df)<-colnames(met)[5:ncol(met)]
num.gene.mat.df.long.1 <- num.gene.mat.df %>% 
  mutate(ThrVal=c(seq(10,100,10))) %>% 
  pivot_longer(1:ncol(num.gene.mat.df), names_to = "Trait", values_to = "GeneNumber")


#for Leaf_Length, Leaf_Width, Leaf Number
#pheno data
pheno <- read.csv("Pheno/21_lettuce_02_pheno.csv", header=T, check.names = FALSE)

num.gene.mat<-matrix(NA,ncol=3,nrow=1)
num.gene.mat.df<-matrix(NA,ncol=3,nrow=0) #or 11
colnames(num.gene.mat) <- colnames(pheno[,c(5,6,8)])
colnames(num.gene.mat.df) <- colnames(pheno[,c(5,6,8)])

#Extract each gene number
for (ThrVal in c(seq(10,100,10))){
  print(ThrVal)
  for (i in 4:6){
    print(paste(ThrVal,i))
    #select phenp
    pheno.sel <- pheno %>% select(1:3,5,6,8)
    
    Trait=names(pheno.sel)[i]
    
    FinalVar<-read.csv(paste(MLtoday,"_",Algo,"/",Trait,"/BorutaAnotation/",ThrVal,"/FinalVar.csv", sep=""), row.names = 1)
    num.gene.mat[1,i-3] <- nrow(FinalVar)
  } 
  num.gene.mat.df <- rbind.data.frame(num.gene.mat.df, num.gene.mat)
}

num.gene.mat.df.long.2 <- num.gene.mat.df %>% 
  mutate(ThrVal=c(seq(10,100,10))) %>% 
  pivot_longer(1:ncol(num.gene.mat.df), names_to = "Trait", values_to = "GeneNumber")

#for Growth
#pheno data
pheno <- read.csv("Pheno/21_lettuce_02_pheno.csv", header=T, check.names = FALSE)

num.gene.mat<-matrix(NA,ncol=1,nrow=1)
num.gene.mat.df<-matrix(NA,ncol=1,nrow=0) #or 11
colnames(num.gene.mat) <- colnames(pheno[4])
colnames(num.gene.mat.df) <- colnames(pheno[4])

#Extract each gene number
for (ThrVal in c(seq(10,80,10))){
  print(ThrVal)
  for (i in 4){
    print(paste(ThrVal,i))
    #select phenp
    pheno.sel <- pheno %>% select(1:4)
    
    Trait=names(pheno.sel)[i]
    
    FinalVar<-read.csv(paste(MLtoday,"_",Algo,"/",Trait,"/BorutaAnotation/",ThrVal,"/FinalVar.csv", sep=""), row.names = 1)
    num.gene.mat[1,i-3] <- nrow(FinalVar)
  } 
  num.gene.mat.df <- rbind.data.frame(num.gene.mat.df, num.gene.mat)
}

num.gene.mat.df.long.3 <- num.gene.mat.df %>% 
  mutate(ThrVal=c(seq(10,80,10))) %>% 
  pivot_longer(1:ncol(num.gene.mat.df), names_to = "Trait", values_to = "GeneNumber")

#for RLW
#pheno data
pheno <- read.csv("Pheno/21_lettuce_02_pheno.csv", header=T, check.names = FALSE)

num.gene.mat<-matrix(NA,ncol=1,nrow=1)
num.gene.mat.df<-matrix(NA,ncol=1,nrow=0) #or 11
colnames(num.gene.mat) <- colnames(pheno[7])
colnames(num.gene.mat.df) <- colnames(pheno[7])

#Extract each gene number
for (ThrVal in c(seq(10,90,10))){
  print(ThrVal)
  for (i in 4){
    print(paste(ThrVal,i))
    #select phenp
    pheno.sel <- pheno %>% select(1:3,7)
    
    Trait=names(pheno.sel)[i]
    
    FinalVar<-read.csv(paste(MLtoday,"_",Algo,"/",Trait,"/BorutaAnotation/",ThrVal,"/FinalVar.csv", sep=""), row.names = 1)
    num.gene.mat[1,i-3] <- nrow(FinalVar)
  } 
  num.gene.mat.df <- rbind.data.frame(num.gene.mat.df, num.gene.mat)
}

num.gene.mat.df.long.4 <- num.gene.mat.df %>% 
  mutate(ThrVal=c(seq(10,90,10))) %>% 
  pivot_longer(1:ncol(num.gene.mat.df), names_to = "Trait", values_to = "GeneNumber")

#Combine
num.gene.mat.df.long <- rbind(num.gene.mat.df.long.1, num.gene.mat.df.long.2, num.gene.mat.df.long.3, num.gene.mat.df.long.4)

#Setting of group turn
num.gene.mat.df.long$Trait <- factor(num.gene.mat.df.long$Trait, levels=c("Chicoric acid", "Chlrogenic acid","Cy3_3MG",
                                                                          "Cy3_6MG","Cy3G","Q3_6MbGG",
                                                                          "Q3_6MGG","Q3_6MG","Q3G",
                                                                          "Growth", "Leaf_Length","Leaf_Number","Leaf_Width","RLW"))

#visualization of gene number
#ggplot
#Fig.6b
ggplot(num.gene.mat.df.long, aes(x=ThrVal, y=GeneNumber))+
  geom_bar(stat="identity", fill="dodgerblue4")+
  geom_text(aes(label=GeneNumber), vjust=-0.2, color="dodgerblue4", size=2.5)+
  #stat_smooth(method="loess",color="red",size=0.5)+
  facet_wrap(~Trait, ncol=7)+
  theme_bw()+
  theme(text=element_text(size=12), 
        axis.text=element_text(size=12, color="black"),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        #axis.title.x = element_blank(),
        )+
  scale_x_continuous(breaks=seq(10,100,10), limits=c(0,110), expand=c(0,0))+
  labs(y="Number of selected genes", x="Number of times selected in 100 repetitions")
ggsave(paste(MLtoday,"_",Algo,"/Fig/",today,"_Barplot_GeneNumber_All.pdf",sep=""), width=10, height=4)


