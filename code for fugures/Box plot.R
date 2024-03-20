library(ggplot2)
library(ggpubr)
library(magrittr)
df <- read.table("111.txt",header = T,sep="\t",na.strings = "NA")
ggplot(df, aes(y=value, x=group,color=group))+facet_grid(~position)+
  geom_boxplot(outlier.shape=NA,color = c("XX","XX","XX"))+
  geom_point(size=1,alpha=0.6,shape=16,position=position_jitter(0.35))+
  theme_bw()+theme(panel.grid=element_blank())+
  scale_color_manual(values = c("XX","XX","XX")) +xlab(NULL)
