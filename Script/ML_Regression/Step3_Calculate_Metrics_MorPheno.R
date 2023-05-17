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

## Modelling based on selected feature----
#Pheno data
pheno <- read.csv("Pheno/Light_MorPheno.csv", header=T, check.names = FALSE)

#for Leaf_Length, Leaf_Width, Leaf Number
#empty matrix for Perform.df
Perform.df <- matrix(NA,nrow=0, ncol=6)
for (i in c(5,6,8)){
  print(i)
  #select phenp
  pheno.sel <- pheno %>% select(1,3,i)
  
  Trait=names(pheno.sel)[3]
  
  for (ThrVal in seq(10,100,10)){
    Performance<-read.csv(paste(today,"_",Algo,"/",Trait,"/Performance/",ThrVal,"/Model_Performance.csv", sep=""), row.names = 1)
    Performance <- Performance %>% 
      mutate(Trait=Trait, ThrVal=ThrVal)
    
    Perform.df <- rbind.data.frame(Perform.df, Performance)
    }
}

#for Growth
i=4
pheno.sel <- pheno %>% select(1,3,i)
  
Trait=names(pheno.sel)[3]
  
  for (ThrVal in seq(10,80,10)){
    Performance<-read.csv(paste(today,"_",Algo,"/",Trait,"/Performance/",ThrVal,"/Model_Performance.csv", sep=""), row.names = 1)
    Performance <- Performance %>% 
      mutate(Trait=Trait, ThrVal=ThrVal)
    
    Perform.df <- rbind.data.frame(Perform.df, Performance)
  }

#for RLW

i=7
pheno.sel <- pheno %>% select(1,3,i)

Trait=names(pheno.sel)[3]

for (ThrVal in seq(10,90,10)){
  Performance<-read.csv(paste(today,"_",Algo,"/",Trait,"/Performance/",ThrVal,"/Model_Performance.csv", sep=""), row.names = 1)
  Performance <- Performance %>% 
    mutate(Trait=Trait, ThrVal=ThrVal)
  
  Perform.df <- rbind.data.frame(Perform.df, Performance)
}

#dir
dir.create(paste(today,"_",Algo,"/PerformanceRawData/Pheno",sep=""))

#output
write.csv(Perform.df, paste(today,"_",Algo,"/PerformanceRawData/Pheno/AllMetrics_Pheno.csv", sep=""))

