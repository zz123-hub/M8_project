library(bruceR)
library(texreg)
df<- read.csv("all_abundance.csv", check.names = F)
df$Gender<- as.factor(df$Gender)
df$Age<- as.factor(df$Age)
X.names<- colnames(df[4:80]) 
M.names<- colnames(df[81:147]) 
Y.names<- colnames(df[148:150]) 
result1<- list()
k = 1
for (g in Y.names){
  for (i in X.names) {
    for (j in M.names) {
      k = 1
      result1[[i]] <- list()
      result1[[i]][[k]] <- PROCESS(df, y = g, x = i,
                                   meds = j, covs = c("Gender","Age"),
                                   ci="boot", nsim=1000, seed =1)
      k = k +1
    }
  }
}
