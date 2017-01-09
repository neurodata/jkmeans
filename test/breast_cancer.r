require("jkmeans")
require("mclust")



setwd("~/git/jkmeans/wdbc/")

wdbc<- read.csv("./wdbc.data.txt",header = F,stringsAsFactors = F)

M<- wdbc$V2
X<- wdbc[,3:ncol(wdbc)]
X<- as.matrix(X)

table(M)
K<- 2

gmm<- jkmeansEM(X,k = K,useKmeansIni = T,meansIni = matrix(rep(1,2*10),nrow = 2),sigma2_ini = 1,fixW = F, tau = 1)

kmeans<- jkmeansEM(X,k = K,useKmeansIni = T,meansIni = matrix(rep(1,2*10),nrow = 2),sigma2_ini = 1,fixW = T, tau = 1000)

jk<- jkmeansEM(X,k = K,useKmeansIni = T,meansIni = matrix(rep(1,2*10),nrow = 2),sigma2_ini = 1,fixW = F,   tau = 1000)

jk_point1<- jkmeansEM(X,k = K,useKmeansIni = T,meansIni = matrix(rep(1,2*10),nrow = 2),sigma2_ini = 1,fixW = F,   tau = 0.1)


gmm$zeta[2,]
kmeans$zeta[2,]
jk$zeta[2,]
jk_point1$zeta[2,]

adjustedRandIndex(gmm$M,M)
adjustedRandIndex(kmeans$M,M)
adjustedRandIndex(jk$M,M)
adjustedRandIndex(jk_point1$M,M)