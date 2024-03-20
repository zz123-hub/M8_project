aaa <- read.table("222.txt",header = T,row.names = 1,sep="\t")
library(ggplot2)
ggplot(aaa,aes(x=PC1,y=PC2,colour=group))+
  geom_point(alpha =.7, size=4)+scale_size_area()+stat_ellipse(level=0.5,mapping= aes(colour=group))+
  scale_colour_brewer(type = "div", palette = "Dark2")+theme_classic()+
  geom_vline(xintercept = 0, color = 'black', size = 0.4, linetype = 4)+
  geom_hline(yintercept = 0, color = 'black', size = 0.4, linetype = 4)+
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title=element_blank(),legend.position = c(0.2,0.9),
        legend.key.size=unit(5,'mm'),legend.key.width=unit(5,'mm'))+
  labs(title="Group A ",x="PC1(15.31%)", y="PC2(4.64%)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("xx","xx"))
