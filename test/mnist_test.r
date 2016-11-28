setwd("~/git/jkmeans/test/")
source("kmpp.r")
require("jkmeans")
require("mclust")
source("PCAC.R")


mnist<- read.csv("../mnist_data/mnist.csv",header = F)

pick<- c(0,3,7,9)
# pick<- c(0:9)

subset<- c(sapply(pick*1000, function(x){x+ c(1:100)*5}))

# subset<- c(sapply(pick*1000, function(x){x+ c(1:1000)}))

mnist_subset<- mnist[,subset]
Y<-t(mnist_subset)


M<- rep(pick,each=100)

# M<- rep(pick,each=1000)
K<- 4

# mclust<-Mclust(Y,G = K,modelNames = "VEI")
# raw_M<-apply(mclust$z, 1, function(x){c(1:K)[x==max(x)]})

# svd0<- svd(Y)

######

#######


mclust<-Mclust(Y,G = K)
raw_M<-apply(mclust$z, 1, function(x){c(1:K)[x==max(x)]})
adjustedRandIndex(raw_M, M)

ari<- numeric()
for(i in 1:10){
  pca_M<-PCAC(Y,K,d = 6)
  ari<- c(ari,adjustedRandIndex(pca_M$M, M))
}

mean(ari)

arc_ari<- numeric()
for(i in 1:10){
  arc_M<- rDARC(Y,d = 6,k = K,meansIni = matrix(1,2),steps = 1E3,sigma2_ini = 0.1,randomStart = F,ver = 1,fixW = F)
  arc_ari<- c(arc_ari,  adjustedRandIndex(arc_M$M, M))
}

image(matrix(arc_M$mu[1,]%*%arc_M$EV,28,28))
image(matrix(arc_M$mu[2,]%*%arc_M$EV,28,28))
image(matrix(arc_M$mu[3,]%*%arc_M$EV,28,28))
image(matrix(arc_M$mu[4,]%*%arc_M$EV,28,28))

plot(arc_M$M)


arc_M$X[4,]
arc_M$mu[4,]

image(matrix(arc_M$X[1,]%*%arc_M$EV,28,28))
image(matrix(arc_M$mu[1,]%*%arc_M$EV,28,28))
image(matrix(Y[1,],28,28))

image(matrix(arc_M$X[101,]%*%arc_M$EV,28,28))
image(matrix(arc_M$mu[2,]%*%arc_M$EV,28,28))
image(matrix(Y[101,],28,28))

image(matrix(arc_M$X[201,]%*%arc_M$EV,28,28))
image(matrix(arc_M$mu[4,]%*%arc_M$EV,28,28))
image(matrix(Y[201,],28,28))

image(matrix(arc_M$X[350,]%*%arc_M$EV,28,28))
image(matrix(arc_M$mu[3,]%*%arc_M$EV,28,28))
image(matrix(Y[303,],28,28))





mean(arc_ari)


plot(raw_M)
plot(pca_M$M)



