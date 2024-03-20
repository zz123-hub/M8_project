library(ggplot2)
aaa <- read.table("111.txt",header = T,row.names = 1,sep="\t")
p <-ggplot(aaa,aes(x=PCoA1,y=PCoA2,colour=tt))+
  geom_point(size=4)+scale_size_area()+stat_ellipse(level=0.5,mapping= aes(colour=group))+
  scale_colour_brewer(type = "div", palette = "Dark2")+theme_classic()+
  geom_vline(xintercept = 0, color = 'black', size = 0.4, linetype = 4)+
  geom_hline(yintercept = 0, color = 'black', size = 0.4, linetype = 4)+
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title=element_blank(),legend.position = c(0.1,0.73),
        legend.key.size=unit(3,'mm'),legend.key.width=unit(3,'mm'))+
  labs(x="PC1(39.44%)", y="PC2(16.08%)")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_colour_manual(values = c("xx","xx"))
p+x11(width=6.5,height=5.5)