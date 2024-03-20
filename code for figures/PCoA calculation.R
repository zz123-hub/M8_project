library(vegan)
library(permute)
library(lattice)
library(ggplot2)
library(plyr)
dis <- read.delim('distance.txt', row.names = 1, sep = '\t', stringsAsFactors = F, check.names = F)
group <- read.delim('group.txt', sep = '\t', stringsAsFactors = F)
pcoa <- cmdscale(as.dist(dis), k = (nrow(dis) - 1), eig = T)
pcoa$eig
point <- data.frame(pcoa$point)
write.csv(point, 'pcoa.sample.csv')
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
pcoa_eig
sample_site <- data.frame({pcoa$point})[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')
sample_site <- merge(sample_site,group, by = 'names', all.x = T)
write.csv(sample_site, 'sample_site.csv',quote = F)
extract_letter <- function(x){
  B<-strsplit(x,"") 
  C<-""
  for (i in 1:length(B[[1]])){
    if (is.na(as.numeric(B[[1]][i]))){
      C <- paste(C,B[[1]][i],sep="")
    }
    else{
      break
    }
  }
  return(C)
}