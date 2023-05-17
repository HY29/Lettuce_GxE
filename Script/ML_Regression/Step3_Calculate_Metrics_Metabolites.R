#########################################################
### Calculation of R and RMSE value in genomic prediction
#########################################################

library(Metrics)
library(tidyverse)

##SetParameters
#Set seed
set.seed(123)
#Current date
#today <- Sys.Date()
today <- "XXXX-XX-XX"
print(today,quote=F)
#Algorithms
Algo=print("RF_Boruta")
#Number of repetition
Repeat=100 
#Threshhold for feature selection in 100 Rep.
ThrVal=100

## Modelling based on selected feature----
#Pheno data
pheno <- read.csv("Pheno/Light_Metabolites.csv", header=T, check.names = FALSE)

R2.mat<-matrix(NA,ncol=ncol(pheno[,5:ncol(pheno)]),nrow=Repeat)
RMSE.mat<-matrix(NA,ncol=ncol(pheno[,5:ncol(pheno)]),nrow=Repeat)
MAE.mat<-matrix(NA,ncol=ncol(pheno[,5:ncol(pheno)]),nrow=Repeat)

RepNo <- paste("Rep_",1:Repeat,"", sep="")

rownames(R2.mat) <- RepNo
colnames(R2.mat) <- colnames(pheno[,5:ncol(pheno)])
rownames(RMSE.mat) <- RepNo
colnames(RMSE.mat) <- colnames(pheno[,5:ncol(pheno)])
rownames(MAE.mat) <- RepNo
colnames(MAE.mat) <- colnames(pheno[,5:ncol(pheno)])

for (ThrVal in c(1, seq(10,100,10))){
  print(ThrVal)

for (i in 5:ncol(pheno)){
  print(i)
  #select phenp
  pheno.sel <- pheno %>% select(1,3,i)
  
  Trait=names(pheno.sel)[3]
  
  for (z in 1:Repeat){
    print(paste(i,Trait,z))
  　Test<-read.csv(paste(today,"_",Algo,"/",Trait,"/Test/TestData_Rep",z,".csv", sep=""), row.names = 1)
    Pred<-read.csv(paste(today,"_",Algo,"/",Trait,"/Prediction/",ThrVal,"/Prediction_Rep",z,".csv", sep=""), row.names = 1)
    plot(Test[,ncol(Test)], Pred[,1])
    
    R2.mat[z,i-3] <- cor(Test[,ncol(Test)], Pred[,1])^2
    RMSE.mat[z,i-3] <- rmse(Test[,ncol(Test)], Pred[,1])
    MAE.mat[z,i-3] <- mae(Test[,ncol(Test)], Pred[,1])
  }
}

#dir
dir.create(paste(today,"_",Algo,"/PerformanceRawData", sep=""))
dir.create(paste(today,"_",Algo,"/PerformanceRawData/",ThrVal,sep=""))

#output
write.csv(R2.mat, paste(today,"_",Algo,"/PerformanceRawData/",ThrVal,"/AllMet_R2.csv", sep=""))
write.csv(RMSE.mat, paste(today,"_",Algo,"/PerformanceRawData/",ThrVal,"/AllMet_RMSE.csv", sep=""))
write.csv(MAE.mat, paste(today,"_",Algo,"/PerformanceRawData/",ThrVal,"/AllMet_MAE.csv", sep=""))
}
