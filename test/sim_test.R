setwd("~/git/jkmeans/test/")
source("kmpp.r")
require("jkmeans")
require("mclust")
source("PCAC.R")


mnist<- read.csv("../mnist_data/mnist.csv",header = F)

# pick<- c(0,3,7,9)
pick<- c(1,3)

subset<- c(sapply(pick*1000, function(x){x+ c(1:100)}))

# subset<- c(sapply(pick*1000, function(x){x+ c(1:1000)}))

mnist_subset<- mnist[,subset]
Y<-t(mnist_subset)


M<- rep(pick,each=100)

# M<- rep(pick,each=1000)
K<- length(pick)

# mclust<-Mclust(Y,G = K,modelNames = "VEI")
# raw_M<-apply(mclust$z, 1, function(x){c(1:K)[x==max(x)]})

# svd0<- svd(Y)

######

#######


# mclust<-Mclust(Y,G = K)
# raw_M<-apply(mclust$z, 1, function(x){c(1:K)[x==max(x)]})
# adjustedRandIndex(raw_M, M)

#ari<- numeric()
#for(i in 1:10){
#  pca_M<-PCAC(Y,K,d = 6)
#  ari<- c(ari,adjustedRandIndex(pca_M$M, M))
#}
#max(ari)

resultPCA<- lapply(c(1:10),function(x) PCAC(Y,K,d = 6))

resultARC<- lapply(c(1:10),function(x) rDARC(Y,d = 6,k = K,meansIni = matrix(1,2),steps = 1E3,sigma2_ini = 0.1,randomStart = F,ver = 1))


l_loglik<- (sapply(resultARC, function(x)x$loglik))
arc_M<- resultARC[[c(1:10)[l_loglik== min(l_loglik)]]]


# save(resultPCA,file="resultPCAMNIST.Rda")
# save(resultARC,file="resultARCMNIST.Rda")

# arc_M$loglik

image(matrix(arc_M$mu[1,]%*%arc_M$EV,28,28))
image(matrix(arc_M$mu[2,]%*%arc_M$EV,28,28))


image(matrix(arc_M$X[1,]%*%arc_M$EV,28,28))
image(matrix(arc_M$X[2,]%*%arc_M$EV,28,28))

plot(arc_M$X,col=arc_M$M+1)

#image(matrix(arc_M$mu[3,]%*%arc_M$EV,28,28))
#image(matrix(arc_M$mu[4,]%*%arc_M$EV,28,28))
#
# 
# image(matrix(arc_M$mu[1,]%*%arc_M$EV,28,28))
# image(matrix(arc_M$mu[2,]%*%arc_M$EV,28,28))
# image(matrix(arc_M$mu[3,]%*%arc_M$EV,28,28))
# image(matrix(arc_M$mu[4,]%*%arc_M$EV,28,28))
# image(matrix(arc_M$mu[5,]%*%arc_M$EV,28,28))
# image(matrix(arc_M$mu[6,]%*%arc_M$EV,28,28))
# image(matrix(arc_M$mu[7,]%*%arc_M$EV,28,28))
# image(matrix(arc_M$mu[8,]%*%arc_M$EV,28,28))
# image(matrix(arc_M$mu[9,]%*%arc_M$EV,28,28))
# image(matrix(arc_M$mu[10,]%*%arc_M$EV,28,28))

#image(matrix(Y[1,],28,28))
#
# image(matrix(arc_M$X[101,]%*%arc_M$EV,28,28))
#image(matrix(arc_M$mu[2,]%*%arc_M$EV,28,28))
#image(matrix(Y[101,],28,28))
#
#image(matrix(arc_M$X[201,]%*%arc_M$EV,28,28))
#image(matrix(arc_M$mu[4,]%*%arc_M$EV,28,28))
#image(matrix(Y[201,],28,28))
#
#image(matrix(arc_M$X[350,]%*%arc_M$EV,28,28))
#image(matrix(arc_M$mu[3,]%*%arc_M$EV,28,28))
#image(matrix(Y[303,],28,28))

# load("resultARCMNIST.Rda")
# load("resultPCAMNIST.Rda")

max(sapply(  resultPCA, function(x) adjustedRandIndex(x$M, M)))
max(sapply(  resultARC, function(x) adjustedRandIndex(x$M, M)))
#
plot(resultARC[[1]]$X,col=resultARC[[1]]$M+1)
plot(resultARC[[2]]$X,col=resultARC[[2]]$M+1)

##
mean(sapply(c(1:10), function(i){
  mat_M <- matrix(resultARC[[i]]$M,nrow=100)
  sum(apply(mat_M, MARGIN = 2, function(x){
    tb<- table(x)
    sum(x != names(tb)[tb==max(tb)])
  }))/1000
}))
##
mean(sapply(c(1:10), function(i){
  mat_M <- matrix(resultPCA[[i]]$M,nrow=100)
  sum(apply(mat_M, MARGIN = 2, function(x){
    tb<- table(x)
    sum(x != names(tb)[tb==max(tb)])
  }))/1000
}))
##

#
load("mclust_MNIST.Rda")

raw_Mat<- matrix(raw_M,nrow=100)
sum(apply(raw_Mat, MARGIN = 2, function(x){
  tb<- table(x)
  sum(x != names(tb)[tb==max(tb)])
}))/1000
#

adjustedRandIndex(raw_M,M)
