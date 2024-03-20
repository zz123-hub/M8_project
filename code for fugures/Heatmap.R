library(pheatmap)
all<-read.table("111.txt",header=TRUE,sep="\t",row.names=1)
group<-as.factor(all[,1])
data<-all[,-1]
annotation<-data.frame(group) 
rownames(annotation)<-rownames(data)
pheatmap(t(data),border=FALSE,
         annotation=annotation,
         cellwidth=10,cellheight=10,cluster_col = FALSE,height = 10,
         cluster_rows = FALSE,
         fontsize_row=12,breaks = seq(-5,5,length.out = 100),
         fontsize_col=12,filename="111.pdf",
         color = colorRampPalette(c("#441E51","#374F6F","#108D86","#55A274","#78AD63","#A8BF4B","#CDCF30","#FDE300"))(100))

