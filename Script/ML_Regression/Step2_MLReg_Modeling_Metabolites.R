#########################################################
### Regression modeling by machine learning algorithms
### Transcriptomic prediction 
### Draft script 210203 HY
### Update 220526 HY
#########################################################

#Fig.7a-d, Fig.S8

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
#today <- Sys.Date()
today <- "XXXX-XX-XX"
print(today,quote=F)
#Algorithms
Algo=print("RF_Boruta")
#Check Algorithms
modelLookup("rf")
#Number of hyperparameter
ParaNum=nrow(modelLookup("rf"))
#Number of repetition
Repeat=100 #Minimum=10, Better=50, Best=100

###Datainput----
#Pheno data
pheno <- read.csv("Pheno/Light_Metabolites.csv", header=T, check.names = FALSE)
#Met_No <- paste("Met_",1:(ncol(pheno)-3), sep="")
#colnames(pheno)[4:12] <- Met_No
#Input Gene annotation Inf
GeneAno<-read.csv("GeneInfo/Lsativa_467_v5.annotation_info.txt", header=T, sep="\t")

## Modelling based on selected feature----
for (z in 5:ncol(pheno)){
  print(z)
  #select phenp
  pheno.sel <- pheno %>% select(1,3,z)
  
  Trait=names(pheno.sel)[3]
  
  for (i in 1:Repeat){
  print("START Modelling")
  print(paste(z,i))
  
  #Reconstruct Train & Test data from final variables by Boruta
  Train <- read.csv(paste(today,"_",Algo,"/",Trait,"/Training/TrainingData_Rep",i,".csv", sep=""), row.names = 1)
  Test <- read.csv(paste(today,"_",Algo,"/",Trait,"/Test/TestData_Rep",i,".csv", sep=""), row.names = 1)
  
  #Feature Selection
  SelVar_num_mean <- read.csv(paste(today,"_",Algo,"/",Trait,"/Boruta/Feature_Decision_Mean.csv", sep=""), row.names=1)

  for (ThrVal in c(1, seq(10,100,10))){
    print(ThrVal)
    
    #mkdir
    dir.create(paste(today,"_",Algo,"/",Trait,"/BorutaAnotation", sep=""))
    dir.create(paste(today,"_",Algo,"/",Trait,"/Prediction/",ThrVal, sep=""))
    dir.create(paste(today,"_",Algo,"/",Trait,"/BorutaAnotation/",ThrVal, sep=""))
    dir.create(paste(today,"_",Algo,"/",Trait,"/Performance/",ThrVal, sep=""))
    dir.create(paste(today,"_",Algo,"/",Trait,"/Importance/",ThrVal, sep=""))
    dir.create(paste(today,"_",Algo,"/",Trait,"/Fig/",ThrVal, sep=""))
    dir.create(paste(today,"_",Algo,"/",Trait,"/PDP/",ThrVal, sep=""))
    
    #FinalVariant selection
    FinalVar <- rownames(SelVar_num_mean %>% filter(Mean>=ThrVal)) #we can flexibilily change filtering condition. 
  
    Train.Boruta <- Train[,FinalVar] %>% bind_cols(Train[,ncol(Train)])
    colnames(Train.Boruta)[ncol(Train.Boruta)] <- Trait
    Test.Boruta <- Test[,FinalVar]  %>% bind_cols(Test[,ncol(Test)])
    colnames(Test.Boruta)[ncol(Test.Boruta)] <- Trait
  
    ##StepXX_Annotation of selected genes
    #Combine FinalVar and Annotation
    FinalVar <- FinalVar %>% as.data.frame()
    colnames(FinalVar)<- "GeneID"
    FinalVar_Ano <- left_join(FinalVar, GeneAno,by="GeneID", na_matches="never")
    FinalVar_Ano_Gene <- FinalVar_Ano %>% distinct(GeneID, .keep_all=TRUE)
    
    #Output
    write.csv(FinalVar, paste(today,"_",Algo,"/",Trait,"/BorutaAnotation/",ThrVal,"/FinalVar.csv", sep=""))
    write.csv(FinalVar_Ano, paste(today,"_",Algo,"/",Trait,"/BorutaAnotation/",ThrVal,"/FinalVar_GeneAnotation.csv", sep=""), row.names = FALSE)
    write.csv(FinalVar_Ano_Gene, paste(today,"_",Algo,"/",Trait,"/BorutaAnotation/",ThrVal,"/FinalVar_GeneAnotation_GeneLevel.csv", sep=""), row.names = FALSE)
    
    #EmptyMatrix for ModelPerform 
    ModelPerform_DF <- matrix(NA,nrow=Repeat,ncol=3)
    colnames(ModelPerform_DF) <- c("RMSE","R2","MAE")
    #EmptyMatrix for ModelImportance
    ModelVarImp_DF <- matrix(NA,nrow=length(FinalVar), ncol=Repeat)
    #EmptyMatrix for BestHyperParameters
    BestParameters_DF <- matrix(NA,nrow=ParaNum, ncol=Repeat)
    #EmptyMatrix for pdp
    pdp_DF <- matrix(NA,nrow=0, ncol=4)
    colnames(pdp_DF) <- c("X","yhat","Variables","Rep")

    ###Step2_Modeling with hyperparameter Tuning by RandomSearch-----
    fitControl <- trainControl(method = "repeatedcv", 
                             number = 10, 
                             repeats = 10, 
                             search="random",
                             verboseIter = TRUE)

    BestModel <- train(y = Train.Boruta[,ncol(Train.Boruta)], 
                     x = Train.Boruta %>% select(-.data[[Trait]]),
                     method = "rf", 
                     trControl = fitControl, 
                     verbose =  TRUE,
                     tuneLength = 30)
    
    BestParameters <- data.frame(BestModel$bestTune)
    rownames(BestParameters_DF) <- colnames(BestParameters)
    BestParameters_DF[1,i] <- BestParameters[,1]

    ###Step3_Prediction by BestModel----
    Pred <- predict(BestModel, newdata=Test.Boruta %>% select(-.data[[Trait]]))
    ModelPerform <- postResample(pred = Pred, obs = Test.Boruta[,ncol(Test.Boruta)])
    plot(Pred, Test.Boruta[,ncol(Test.Boruta)])
    
    write.csv(Pred, paste(today,"_",Algo,"/",Trait,"/Prediction/",ThrVal,"/Prediction_Rep",i,".csv", sep=""))
  
    ModelPerform_DF[i,1] <- ModelPerform[1] #RMSE
    ModelPerform_DF[i,2] <- ModelPerform[2] #R2
    ModelPerform_DF[i,3] <- ModelPerform[3] #MAE
  
    ###Step4_Variable Importance----
    Imp <- varImp(BestModel, useModel = TRUE, scale=FALSE)
    rownames(ModelVarImp_DF) <- rownames(Imp$importance)
    ModelVarImp_DF[,i] <- Imp$importance[rownames(ModelVarImp_DF),]
  
    ###Step5_Partial Dependence Plot (PDP)----
    for (j in 1:(ncol(Train.Boruta)-1)){
    pdp <- pdp::partial(BestModel, 
                      train=Train.Boruta, 
                      pred.var=paste(colnames(Train.Boruta)[j]),
                      grid.resolution=NULL, #default value=51
                      type="regression") %>%
    as.data.frame()
    
    pdp <- cbind.data.frame(pdp,rep(colnames(Train.Boruta)[j], nrow(pdp)),paste("Rep_",i,"",sep=""))
    colnames(pdp) <- c("X","yhat","Variables","Rep")
    pdp_DF <- rbind.data.frame(pdp_DF,pdp)
  }

###Step6_Output of ModelPerformance----
#Output Raw ModelPerformance
RepNo <- paste("Rep_",1:Repeat,"", sep="")
ModelPerform_DF <- cbind.data.frame(RepNo,ModelPerform_DF)
colnames(BestParameters_DF) <- RepNo
colnames(ModelVarImp_DF) <- RepNo
write.csv(BestParameters_DF, paste(today,"_",Algo,"/",Trait,"/Performance/",ThrVal,"/BestParameters.csv", sep=""))
write.csv(ModelPerform_DF, paste(today,"_",Algo,"/",Trait,"/Performance/",ThrVal,"/Model_Performance.csv", sep=""))
write.csv(ModelVarImp_DF, paste(today,"_",Algo,"/",Trait,"/Importance/",ThrVal,"/Model_Impotance.csv", sep=""))
#Error_ntree_Plot
pdf(paste(today,"_",Algo,"/",Trait,"/Fig/",ThrVal,"/BestModel_Plot.pdf", sep=""), width=4, height=4)
plot(BestModel$finalModel)
dev.off()

#Output Mean&SD for ModelPerformance
#ModelPerformance
MeanPer <- cbind(mean(ModelPerform_DF[,2]), mean(ModelPerform_DF[,3]), mean(ModelPerform_DF[,4]))
SDPer <- cbind(sd(ModelPerform_DF[,2]),sd(ModelPerform_DF[,3]),sd(ModelPerform_DF[,4]))
PerMean_DF <- rbind(MeanPer, SDPer)
colnames(PerMean_DF) <- c("RMSE","R2","MAE")
rownames(PerMean_DF) <- c("Mean","SD")
write.csv(PerMean_DF,paste(today,"_",Algo,"/",Trait,"/Performance/",ThrVal,"/Model_Performance_MeanSD.csv", sep=""))
#ModelImportance
ImpMean_DF <- matrix(NA,nrow=length(Train.Boruta[,-1]), ncol=2)
colnames(ImpMean_DF) <- c("Mean","SD")
rownames(ImpMean_DF) <- rownames(Imp$importance)
for (i in 1:length(Train.Boruta[,-1])){
  ImpMean_DF[i,1]<- mean(ModelVarImp_DF[i,])
  ImpMean_DF[i,2]<- sd(ModelVarImp_DF[i,])
}
ImpMean_DF <- ImpMean_DF[order(ImpMean_DF[,1], decreasing = T),]
write.csv(ImpMean_DF, paste(today,"_",Algo,"/",Trait,"/Importance/",ThrVal,"/Model_Importance_MeanSD.csv", sep=""))  
#partial dependence plot
write.csv(pdp_DF, paste(today,"_",Algo,"/",Trait,"/PDP/",ThrVal,"/PartialDependencePlot_RawValues.csv", sep=""),row.names = FALSE)  

###Step7_Visualization----
#Extraction Best&Worst model
RepNo <- 1:Repeat
ModelPerform_DF <- cbind.data.frame(RepNo,ModelPerform_DF)
Best <- ModelPerform_DF[ModelPerform_DF$R2==max(ModelPerform_DF$R2),][,1]
Worst <- ModelPerform_DF[ModelPerform_DF$R2==min(ModelPerform_DF$R2),][,1]

#BestModel_Plot
#Data arrangement
BestPlot <- cbind.data.frame(
  read.csv(paste(today,"_",Algo,"/",Trait,"/Test/TestData_Rep",print(Best),".csv", sep=""), row.names=1)[,ncol(Test)],
  read.csv(paste(today,"_",Algo,"/",Trait,"/Prediction/Prediction_Rep",print(Best),".csv", sep=""), row.names=1)
)
colnames(BestPlot) <- c("Observed","Predicted")
#ggplot2
ggplot(BestPlot, aes(y=Predicted, x=Observed))+
  guides(color=FALSE, fill=FALSE)+
  #geom_point(size=3.5, shape=21, fill="blue", color="black", alpha=0.8)+
  geom_abline(slope=1,intecept=0, color="red",size=0.5)+
  geom_point(size=2.0, shape=1, color="blue")+
  #stat_smooth(method="lm", se = T, size =0.6,color="red")+
  theme_bw()+
  theme(text=element_text(size=10),axis.text=element_text(size=10,color="black"),title=element_text(size=8))+
  #theme_cowplot(font_size = 16, line_size = 0.5)+
  labs(title=paste("Best_",Algo,sep=""), 
       subtitle=paste("R2=",format(max(ModelPerform_DF$R2), digits=3),sep=""),
       y="Predicted value",x="Observed value")

ggsave(paste(today,"_",Algo,"/",Trait,"/Fig/",ThrVal,"/BestPlot.pdf", sep=""), width=2, height=2.5)

#WorstModel_Plot
#Data arrangement
WorstPlot <- cbind.data.frame(
  read.csv(paste(today,"_",Algo,"/",Trait,"/Test/TestData_Rep",print(Worst),".csv", sep=""), row.names=1)[,ncol(Test)],
  read.csv(paste(today,"_",Algo,"/",Trait,"/Prediction/Prediction_Rep",print(Worst),".csv", sep=""), row.names=1)
)
colnames(WorstPlot) <- c("Observed","Predicted")
#ggplot2
ggplot(WorstPlot, aes(y=Predicted, x=Observed))+
  guides(color=FALSE, fill=FALSE)+
  #geom_point(size=3.5, shape=21, fill="blue", color="black", alpha=0.8)+
  geom_abline(slope=1,intecept=0, color="red",size=0.5)+
  geom_point(size=2.0, shape=1, color="blue")+
  #stat_smooth(method="lm", se = T, size =0.6,color="red")+
  theme_bw()+
  theme(text=element_text(size=10),axis.text=element_text(size=10,color="black"),title=element_text(size=8))+
  #theme_cowplot(font_size = 16, line_size = 0.5)+
  labs(title=paste("Worst_",Algo,sep=""), 
       subtitle=paste("R2=",format(min(ModelPerform_DF$R2), digits=3),sep=""),
       y="Predicted value",x="Observed value")

ggsave(paste(today,"_",Algo,"/",Trait,"/Fig/",ThrVal,"/WorstPlot.pdf", sep=""), width=2, height=2.5)

#Importance_BarPlot
ImpPlot <- read.csv(paste(today,"_",Algo,"/",Trait,"/Importance/",ThrVal,"/Model_Importance_MeanSD.csv", sep=""))
colnames(ImpPlot) <- c("Variables","Mean","SD")
ImpPlot$Variables <- factor(ImpPlot$Variables, levels=ImpPlot$Variables)
#ggplot2
ggplot(ImpPlot, aes(y=Mean, x=Variables))+
  guides(color=FALSE, fill=FALSE)+
  geom_bar(stat="identity",fill="cyan4",color="black", width=0.7,position=position_dodge(0.6))+
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD, width=5),position=position_dodge(0.6),size=0.5, width=0.4)+
  theme_bw()+
  theme(text=element_text(size=10),axis.text=element_text(size=10,color="black"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  #theme_cowplot(font_size=16, line_size=0.5)+
  #scale_y_continuous(breaks=seq(0, 20, 5), limits=c(0,20), expand = c(0,0))+
  #theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
  labs(title=paste("",Algo,sep=""),y="Importance", x="Features")

ggsave(paste(today,"_",Algo,"/",Trait,"/Fig/",ThrVal,"/ImportanceBarPlot.pdf", sep=""), width=3, height=3)

#partial dependence plot
ggplot(pdp_DF,aes(x=X,y=yhat,color=Rep, group=Rep))+
  geom_line(aes(color=Rep,group=Rep),color="red",size=0.3, alpha=0.8)+
  #stat_summary(fun=mean,geom="line",size=0.1,color="blue")+
  #geom_line(stat="identity",color="blue",size=0.5)+
  #geom_bar(stat="identity",fill="cyan4",color="black", width=0.7,position=position_dodge(0.6))+
  #geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD, width=5),position=position_dodge(0.6),size=0.5, width=0.4)+
  stat_smooth(aes(group=1),method="loess", color="blue", fill="blue",size=0.6)+
  facet_wrap(~Variables, ncol=5, scales="free")+
  theme_bw()+
  theme(text=element_text(size=10),axis.text=element_text(size=10,color="black"),
        strip.text.x = element_text(size=10))+
  #scale_x_continuous(breaks=seq(-2,2,1),limits=c(-2,2),expand=c(0,0))+
  labs(y="Partial dependence",x="Feature values")
ggsave(paste(today,"_",Algo,"/",Trait,"/Fig/",ThrVal,"/PartialDependencePlot_All.pdf", sep=""), width=10, height=15)

  }
  }
  }
