library(ggord)
library(ggplot2)
pca_data=read.table("111.txt",header=T,sep="\t",row=1)
pca1=prcomp(pca_data,center=TRUE,retx=T,scale=T)
pca1
summary(pca1)
plot(pca1$x, main="after PCA")
head(pca1$x)
write.csv(pca1$x,'pc.csv', quote = FALSE)
rm1<-cor(pca_data)
rs1<-eigen(rm1)
val <- rs1$values
(Standard_deviation <- sqrt(val))
(Proportion_of_Variance <- val/sum(val)) 

