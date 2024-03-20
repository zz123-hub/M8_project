con <- file("111.txt","r")
line<-readLines(con,n=1)
line<-unlist(strsplit(line,"\t"))
result<-data.frame()
number<-length(line)
sink("123.txt",append=TRUE)
while( length(line) != 0 ) {
  line<-unlist(strsplit(line,"\t"))
  result[1,1]<-line[1]
  result[1,2]<-wilcox.test(as.numeric(line[2:44]),as.numeric(line[45:85]))$p.value
  write.table(result,sep="\t",quote=FALSE,row.names=F,col.names=F)
  line=readLines(con,n=1)
}
close(con)
sink()
