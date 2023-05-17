##Fig.2b

#Boxplot by ggplot2

#library
library(tidyverse)

#data import
d <- read.csv("Pheno/Light_Metabolites.csv",header=TRUE, check.names=FALSE)

d.long <- d %>% 
  filter(Color=="Red") %>% 
  pivot_longer(5:(length(d)), names_to="Pheno", values_to="values") 


#Setting of group turn
d.long$Light<- factor(d.long$Light, levels=c("WhiteLED","RedLED","BlueLED","FL"))
d.long$Genotype<- factor(d.long$Genotype, levels=c("RLL4","rll4"))
d.long$Cultivar<- factor(d.long$Cultivar, levels=c("RButt","RLeaf","ROak","RCos","FancR","HandR","RFire",
                                         "SunB","WineD","GOak","HandG","GLeaf","GButt","FrIce"))

#ggplot2 
p <- ggplot(d.long, aes(y=values, x=Light, fill=Genotype))+
  #guides(color=FALSE, fill=FALSE, shape=FALSE)+
  stat_boxplot(color="black", geom="errorbar", width=0.6, position = position_dodge(width=0.6))+
  geom_boxplot(color="black", width=0.5,position = position_dodge(width=0.6), outlier.colour = NA)+
  geom_jitter(color="black", shape=21, size=0.8, position = position_dodge(width=0.6))+
  facet_wrap(~ Pheno, ncol=3, scales="free")+
  theme_bw()+
  theme(text = element_text(size = 10), axis.text = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        #axis.text.x = element_text(size=12, color="black", angle=90, hjust=1),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 11, face = "bold"),
        strip.text.y = element_text(angle=0),
        legend.title = element_text(size = 8, face = "plain"),
        legend.text = element_text(size = 8, face = "italic"))+
  scale_color_brewer(palette="Set2")+
  scale_fill_brewer(palette="Set2")+
  labs(y="mg/100 g FW", x="Light conditions")

#plot
print(p)

ggsave("Fig/Fig2b_Light_Metabolites_RLL4.pdf", width=6.0, height=4)

