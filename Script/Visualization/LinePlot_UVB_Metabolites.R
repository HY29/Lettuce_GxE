##Fig.3a

#Line graph by ggplot2

#library
library(tidyverse)
library(plotly)

#data import
d <- read.csv("Pheno/UVB_Metabolites.csv",header=TRUE, check.names=FALSE)

d.long <- d %>% 
  pivot_longer(5:(length(d)), names_to="metabolites", values_to="values") %>% 
  group_by(Light, Cultivar, Color, metabolites) %>% 
  summarize(mean=mean(values), sd=sd(values))


#Setting of group turn
d.long$Light<- factor(d.long$Light, levels=c("Control","UVB"))
d.long$Cultivar<- factor(d.long$Cultivar, levels=c("RLeaf","ROak","FancR","HandR","WineD",
                                         "GLeaf","GButt","FrIce"))


#ggplot2 
p <- ggplot(d.long, aes(y=mean, x=Light, color=Color, fill=Color, shape=Cultivar, group=Cultivar))+
  #guides(color=FALSE, fill=FALSE, shape=FALSE)+
  geom_line(size=0.5)+
  geom_point(size=1.5)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position=position_dodge(0),size=0.5, width=1, alpha=0.8)+
  facet_wrap(~ metabolites, ncol=3, scales="free")+
  theme_bw()+
  theme(text = element_text(size = 10), axis.text = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        #axis.text.x=element_text(size=12, color="black", angle=90, hjust=1),
        axis.text.x = element_blank(),
        #strip.background = element_blank(),
        #strip.text = element_text(size = 11, face = "bold"),
        strip.text.y = element_text(angle=0))+
  scale_shape_manual(values=c(0,1,2,3,4,5,6,7,8,9,15,16,17,18))+
  scale_color_manual(values=c("green","deeppink2"))+
  labs(y="mg/100 g FW", x="Light conditions")

#ggplontly  
ggplotly(p)

#plot
print(p)

ggsave("Fig/Fig3a_UVB_Metabolites.pdf", width=5.5, height=4)

