require("jkmeans")
require("mclust")
setwd("~/git/jkmeans/test/")
source("kmpp.r")
source("PCAC.R")


Y<- matrix(runif(10000,-10,10),ncol=2)

ring<-function(a,b){
  x1<- Y[,1]
  x2<- Y[,2]
  (x1^2+x2^2)>a^2 & (x1^2+x2^2)<b^2
}

M<- ring(0,1)*1 + ring(2,4)*2 + ring(5,7)*3
subset<- (ring(0,1)|ring(2,4)|ring(5,7))

Y<- Y[  subset,]
M<- M[subset]
plot(Y,col=M)

###################################

K<- 3

###############################




mclust<-Mclust(Y,G = K,modelNames = "VEI")
raw_M<-apply(mclust$z, 1, function(x){c(1:K)[x==max(x)]})


# pca_M<-PCAC(Y,K,d=2)

arc_M<- rDARC(Y,d = 1,k = K,meansIni = matrix(1,2),steps = 1E3,sigma2_ini = 0.1,randomStart = T,ver = 1)
# arc_M<- ARC(Y,K,d = 2)

testARC<- lapply(c(1:10),function(x) rDARC(Y,d = 2,k = K,meansIni = matrix(1,2),steps = 1E3,sigma2_ini = 0.1,randomStart = F,ver = 1))

l_loglik<- (sapply(testARC, function(x)x$loglik))

arc_M<- testARC[[c(1:10)[l_loglik== min(l_loglik)]]]

arc_M$loglik

pdf("simX.pdf",6,6)
plot(X,col=M,xlab="X1",ylab = "X2",lwd=2)
dev.off()

pdf("simPCAX.pdf",6,6)
plot(pca_M$X,col=pca_M$M+1,xlab="X1",ylab = "X2",lwd=2)
dev.off()

pdf("simARCX.pdf",6,6)
plot(arc_M$X,col=arc_M$M+1,xlab="X1",ylab = "X2",lwd=2,ylim=range(pca_M$X[,2]))
dev.off()

# plot(arc_M$mu,type = "p",lwd=10)

adjustedRandIndex(raw_M,M)
# adjustedRandIndex(pca_M$M, M)
adjustedRandIndex(arc_M$M, M)

plot(Y, col=arc_M$M+1)
plot(Y, col=raw_M+1)
