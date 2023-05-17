#########################################################
### Regression modeling by machine learning algorithms
### Transcriptomic prediction 
#########################################################
#--------------
#Random Forest
#--------------
#library
library(caret)
library(randomForest)
library(Boruta)
library(pdp)
library(vip)
library(mlbench)
library(tidyverse)
library(tidymodels)
library(edgeR)

##SetParameters
#Set seed
set.seed(123)
#Current date
today <- Sys.Date()
print(today,quote=F)
#Algorithms
Algo=print("RF_Boruta")
#Check Algorithms
modelLookup("rf")
#Number of hyperparameter
ParaNum=nrow(modelLookup("rf"))
#Number of repetition
Repeat=100 #Minimum=10, Better=50, Best=100
#Make Dir
dir.create(paste(today,"_",Algo, sep=""))

###Datainput----
#Pheno data
pheno <- read.csv("Pheno/Light_Metabolites.csv", header=T, check.names = FALSE)
#Met_No <- paste("Met_",1:(ncol(pheno)-3), sep="")
#colnames(pheno)[4:12] <- Met_No

for (z in 5:ncol(pheno)){
  print(z)
#select phenp
pheno.sel <- pheno %>% select(1,3,z)

Trait=names(pheno.sel)[3]

#Make Dir
dir.create(paste(today,"_",Algo,"/",Trait, sep=""))
dir.create(paste(today,"_",Algo,"/",Trait,"/TMMnorm", sep=""))
dir.create(paste(today,"_",Algo,"/",Trait,"/Training", sep=""))
dir.create(paste(today,"_",Algo,"/",Trait,"/Test", sep=""))
dir.create(paste(today,"_",Algo,"/",Trait,"/Prediction", sep=""))
dir.create(paste(today,"_",Algo,"/",Trait,"/Performance", sep=""))
dir.create(paste(today,"_",Algo,"/",Trait,"/Importance", sep=""))
dir.create(paste(today,"_",Algo,"/",Trait,"/PDP", sep=""))
dir.create(paste(today,"_",Algo,"/",Trait,"/Boruta", sep=""))
dir.create(paste(today,"_",Algo,"/",Trait,"/Fig", sep=""))

#Gene expression data
count <- read.csv("GeneExp/GeneExpressionMatrix.count.csv", row.names=1, header=T, sep="")

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

#output TMMnorm
write.csv(TMMnorm.log, paste(today,"_",Algo,"/",Trait,"/TMMnorm/TMMnormalized_log2_count.csv", sep=""))

#Matrix for modeling
d <- cbind.data.frame(pheno.sel[,2:3],t(TMMnorm.log))

###Feature selection by Boruta with repetition----
#EmptyMatrix for selected feature by Boruta
SelVar_DF <- matrix(NA,nrow=nrow(TMMnorm.log), ncol=Repeat) %>% as.data.frame()

for (i in 1:Repeat){
  print("START Feature Selection")
  print(paste(z,i))
  
  ###Statified Random Sampling-----
  #Train (70%), Test (30%)
  #SetSeed
  set.seed(0)
  TrainIndex <- createDataPartition(d[,1], p = .7, 
                                  list = FALSE, 
                                  times = Repeat)
  #Train_scaling
  Train_x_raw <- d[TrainIndex[,i],] %>% select(-1,-2)
  Train_x_mean <- apply(Train_x_raw, 2, mean)
  Train_x_sd <- apply(Train_x_raw, 2, sd)
  Train_x <- scale(Train_x_raw) %>% as.data.frame()
  Train_y <- d[TrainIndex[,i],] %>% select(.data[[Trait]]) #objective variable
  Train <- bind_cols(Train_x, Train_y)
  
  #Test_scaling based on Train mean & sd
  Test_x_raw <- d[-TrainIndex[,i],] %>% select(-1,-2)
  Test_x <- data.frame(matrix(NA, nrow = nrow(Test_x_raw), ncol = ncol(Test_x_raw)))
  colnames(Test_x) <- colnames(Test_x_raw)
  for (m in 1:length(Test_x_raw)) {
    Test_x[m] <- (Test_x_raw[m]-Train_x_mean[m])/Train_x_sd[m]
  }
  Test_y <- d[-TrainIndex[,i],] %>% select(.data[[Trait]]) #objective variable
  Test <- bind_cols(Test_x, Test_y)

  write.csv(Train, paste(today,"_",Algo,"/",Trait,"/Training/TrainingData_Rep",i,".csv", sep=""))
  write.csv(Test, paste(today,"_",Algo,"/",Trait,"/Test/TestData_Rep",i,".csv", sep=""))
  
  ###Feature selection by Boruta----
  #Boruta
  Boruta <- Boruta(Train %>% select(-.data[[Trait]]), 
                   Train_y[,1], 
                   getImp = getImpRfZ, 
                   mcAdj = TRUE,
                   maxRuns = 1000,
                   pValue = 0.01,
                   doTrace = 2)
  
  #Summary of Boruta results
  Boruta.Stats <- attStats(Boruta)
  
  Boruta.Stats %>% 
    count(decision)
  Boruta.Stats %>% 
    filter(decision=="Confirmed")
  
  write.csv(Boruta.Stats, paste(today,"_",Algo,"/",Trait,"/Boruta/BorutaStats_Rep",i,".csv", sep=""))
  
  #Extract Borita desicion
  rownames(SelVar_DF) <- colnames(Train %>% select(-.data[[Trait]]))
  SelVar_DF[,i] <- Boruta.Stats %>% select(decision) 
  
  SelVar_num <- SelVar_DF %>% 
    mutate_all(~gsub(., pattern="Rejected", replacement="0")) %>% 
    mutate_all(~gsub(., pattern="Tentative", replacement="0")) %>% 
    mutate_all(~gsub(., pattern="Confirmed", replacement="1")) 
  
}

##Extract selected feature in all Boruta repetition
SelVar_num <- lapply(SelVar_num[1:Repeat],as.numeric) %>% as.data.frame()
SelVar_num_mean <- apply(SelVar_num,1,sum) %>% as.data.frame()
colnames(SelVar_num_mean) <- c("Mean")
rownames(SelVar_num_mean) <- colnames(Train %>% select(-.data[[Trait]]))
FinalVar <- rownames(SelVar_num_mean %>% filter(Mean==Repeat)) #we can flexibilily change filtering condition. 
#FinalVar <- rownames(SelVar_num_mean %>% filter(Mean>0)) #we can flexibilily change filtering condition. 

#output
RepNo <- paste("Rep_",1:Repeat,"", sep="")
colnames(SelVar_DF) <- RepNo
write.csv(SelVar_DF, paste(today,"_",Algo,"/",Trait,"/Boruta/Feature_Decision.csv", sep=""))
write.csv(SelVar_num_mean, paste(today,"_",Algo,"/",Trait,"/Boruta/Feature_Decision_Mean.csv", sep=""))

}


