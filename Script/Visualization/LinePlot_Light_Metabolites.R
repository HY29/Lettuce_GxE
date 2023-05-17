##Fig.2a

#Line graph by ggplot2

#library
library(tidyverse)

#data
d <- read.csv("Pheno/Light_Metabolites.csv",header=TRUE, check.names=FALSE)

d.long <- d %>% 
  pivot_longer(5:(length(d)), names_to="metabolites", values_to="values") %>% 
  group_by(Light, Cultivar, Color, metabolites) %>% 
  summarize(mean=mean(values), sd=sd(values))

#Setting of group turn
d.long$Light<- factor(d.long$Light, levels=c("WhiteLED","RedLED","BlueLED","FL"))
d.long$Cultivar<- factor(d.long$Cultivar, levels=c("RButt","RLeaf","ROak","RCos","FancR","HandR","RFire",
                                         "SunB","WineD","GOak","HandG","GLeaf","GButt","FrIce"))


#ggplot2 
ggplot(d.long, aes(y=mean, x=Light, color=Color, fill=Color, shape=Cultivar, group=Cultivar))+
  guides(color=FALSE, fill=FALSE, shape=FALSE)+
  geom_line(size=0.5)+
  geom_point(size=1.5)+
  facet_wrap(~ metabolites, ncol=3, scales="free")+
  theme_bw()+
  theme(text = element_text(size = 10), axis.text = element_text(size=10, color="black"),
        axis.title.x = element_blank(),
        #axis.text.x=element_text(size=12, color="black", angle=90, hjust=1),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 11, face = "bold"),
        strip.text.y = element_text(angle=0))+
  scale_shape_manual(values=c(0,1,2,3,4,5,6,7,8,9,15,16,17,18))+
  scale_color_manual(values=c("green","deeppink2"))+
  labs(y="mg/100 g FW", x="Light conditions")
  
#ggplontly  
ggplotly(p)

#plot
print(p)

ggsave("Fig/Fig2a_Light_Metabolites.pdf", width=6.0, height=4)
