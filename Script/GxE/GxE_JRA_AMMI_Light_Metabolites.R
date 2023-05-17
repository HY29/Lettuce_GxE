#########################
## GxE analysis
## JRA & AMMI
#########################

#Fig.4, Fig.S6,7

#Current date
today <- Sys.Date()
print(today, quote=F)

#load packages
library(tidyverse)
library(ggrepel)

#data import
d <- read.csv("Pheno/Light_Metabolites.csv",header=TRUE, check.names=FALSE)
d$Cultivar<- factor(d$Cultivar, levels=c("RButt","RLeaf","ROak","RCos","FancR","HandR","RFire",
                                                   "SunB","WineD","GOak","HandG","GLeaf","GButt","FrIce"))
for (j in 5:ncol(d)){

  print(j)

  Trait=names(d)[j]
  
  #make dir
  dir.create(paste("Fig/GxE/",Trait,sep=""))
  
  #ANOVA test for GXE detection
  twANOVA <- aov(d[,j]~Cultivar*Light, data=d)
  summary(twANOVA)
  
  #Interaction plot by ggplot2
  ggplot(d, aes(y=d[,j], x=Cultivar, color=Light, group=Light))+
  #guides(color=FALSE, fill=FALSE)+
  stat_summary(fun.y= mean, geom="line", size=0.5, width=0.8)+
  geom_point(shape=21, size=1.5, position=position_dodge(0.8))+
  theme_bw()+
  theme(text = element_text(size = 10), axis.text = element_text(size=10, color="black"), 
        axis.text.x = element_text(angle = 60, hjust = 1), 
        axis.title.x = element_blank(),
        plot.title = element_text(size = 7))+
  labs(title = Trait, y="Phenotype", x="genotype", color="Environment")
  
  ggsave(paste("Fig/GxE/",Trait,"/",print(today, quate=F),"_InteractionPlot.pdf", sep=""), width=5, height=3)
  
  #data arrange
  d.mean <- d %>% 
  group_by(Light, Cultivar) %>% 
  summarise(mean=mean(.data[[Trait]]))%>% 
  pivot_wider(names_from=Light, values_from=mean) %>% 
  as.data.frame()
  
  #d.mean <- as.matrix(d.mean)
  rownames(d.mean) <- d.mean[,1]
  d.mean <- d.mean[,2:5]
  
  print(d.mean)
  
  #overall mean
  om <- mean(as.matrix(d.mean))
  #genotypic mean
  gm <- apply(d.mean, 1, mean)
  #environmental mean
  em <- apply(d.mean, 2, mean)
  
  
  ####Joing Regression Analysis(JRA)----
  #prepare empty matirix
  g <- rep(NA, nrow(d.mean))
  b <- rep(NA, nrow(d.mean))
  rv <- rep(NA, nrow(d.mean))
  #tj <- em - om # just in case
  
  #linear modeling
  for(i in 1:nrow(d.mean)) {
    print(i)
    model <- lm(as.numeric(d.mean[i, ]) ~ em)
    g[i] <- coef(model)[1]
    b[i] <- coef(model)[2]
    rv[i] <- summary(model)$sigma^2
    }
  
  #Trait values vs. Environmental index(em)
  #ggplot ver. (in progress)
  d.mean.t <- t(d.mean) %>% 
    as.data.frame()
  
  d.mean.t <- bind_cols(em, d.mean.t) 
  colnames(d.mean.t)[1] <- "em"
  
  d.mean.t.long <- d.mean.t %>% 
    pivot_longer(2:(length(d.mean.t)), names_to="metabolites", values_to="values")
  
  col <- rep("blue3", length(b))
  col[b < quantile(b, 3/4)] <- "green3"
  col[b <= quantile(b, 1/4)] <- "red2"
  
  P <- ggplot()+
    geom_point(data=d.mean.t.long, aes(x=em, y=values), shape=1, color="black", size=2.0)+
    theme_bw()+
    theme(text = element_text(size = 14), axis.text = element_text(size=14, color="black"),
          plot.title = element_text(size = 14))+
    labs(title=Trait, y="Traits values", x="Environmental index")
  
  for (i in 1:14){
    P <- P +geom_abline(intercept = g[i], slope=b[i], color=col[i], size=0.5)
  }
  P <- P+geom_abline(intercept=om, slope=1.0, linetype = "dotted",color="black")
  print(P)
  
  
  ggsave(paste("Fig/GxE/",Trait,"/",print(today, quate=F),"_JRA_Regression.pdf", sep=""), width=3, height=3)
  
  #slope (stability) vs. intercept (trait)
  #ggplot ver
  Type <- c(rep("Red",9), rep("Green",5))
  stab.df <- bind_cols(g,b,rownames(d.mean), Type)
  colnames(stab.df) <- c("g","b","Cultivar","Type")
  
  ggplot(stab.df, aes(x=g, y=b, fill=Type))+
    guides(color=FALSE, fill=FALSE)+
    geom_point(shape=21, size=3.0, color="black")+
    geom_label_repel(aes(label=Cultivar, fill=Type), size=3.0, alpha=0.8, color="black")+
    scale_fill_manual(values=c("green","deeppink2"))+
    geom_hline(yintercept = 1.0, linetype = "dotted")+
    theme_bw()+
    theme(plot.title = element_text(size = 14), text = element_text(size = 14), axis.text = element_text(size=14, color="black"))+
    labs(title = Trait, y="Stability (b)", x="Intercept")
  
  ggsave(paste("Fig/GxE/",Trait,"/",print(today, quate=F),"_Plot_Stability_vs_Intercept.pdf", sep=""), width=3, height=3)
  
  ###Additive Main effect Multiplicative Interaction(AMMI) model----
  #matrix including  GXE effect
  ge <- sweep(sweep(d.mean, 1, gm), 2, em) + om
  head(ge)
  
  sge <- svd(ge)
  
  zg <- sge$u %*% diag(sqrt(sge$d))
  ze <- sge$v %*% diag(sqrt(sge$d))
  rownames(zg) <- names(gm)
  rownames(ze) <- names(em)
  colnames(zg) <- colnames(ze) <- paste0("IPCA", 1:length(sge$d))
  
  #ggplot ver.
  ze.df <- bind_cols(as.data.frame(ze),rownames(ze))
  colnames(ze.df)[5] <- c("Env")
  
  zg.df <- bind_cols(as.data.frame(zg),rownames(zg),Type)
  colnames(zg.df)[5:6] <- c("Cultivar","Type")
  
  ggplot()+
    guides(color=FALSE, fill=FALSE)+
    geom_point(data=zg.df, aes(x=IPCA1, y=IPCA2, color=Type), shape=1, size=3.0)+
    scale_color_manual(values=c("green","deeppink2"))+
    scale_fill_manual(values=c("green","deeppink2"))+
    geom_segment(data=ze.df, aes(x=0, y=0, xend=IPCA1, yend=IPCA2),arrow=arrow(length=unit(0.3,"cm"), type = "open"), color="blue", alpha=0.8, size=0.5)+
    geom_label_repel(data=ze.df, aes(x=IPCA1, y=IPCA2, label=Env), size=3, color="blue", alpha=0.8)+
    geom_text_repel(data=zg.df, aes(x=IPCA1, y=IPCA2, label=Cultivar, color=Type), size=3.0, alpha=0.8, max.overlaps = Inf)+
    theme_bw()+
    theme(text = element_text(size = 14), axis.text = element_text(size=14, color="black"),
          plot.title = element_text(size = 14))+
    geom_hline(yintercept = 0, linetype = "dotted")+
    geom_vline(xintercept = 0, linetype = "dotted")+
    labs(title = Trait, y="IPCA2", x="IPCA1")
  
  ggsave(paste("Fig/GxE/",Trait,"/",print(today, quate=F),"_IPCPlot_AMMI.pdf", sep=""), width=3, height=3)

}

