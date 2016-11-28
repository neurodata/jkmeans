require("jkmeans")
require("mclust")
setwd("~/git/jkmeans/test/")
source("kmpp.r")
source("PCAC.R")


Y<- matrix(rnorm(1000),nrow = 100)

mu<- rep(c(1:3),each=100)

Y<- matrix(rnorm(300*3,mu, sd = 0.5),300)

# fit<- rDARC(Y,d = 2,k = 2,meansIni = matrix(1,2))

# kmpp(Y,3)


#


n<- 300
p<- 1000
d<- 2
sigma2<- 1
tau<- 0.1
Sigma<- diag(tau,d)

###################################

K<- 3
mu1<- matrix(K,n,d)
# randomMu<- unlist(sapply(1:K, function(x){rep(x, runif(1,n/10,n))}))
# mu1[1:length(randomMu)]<- randomMu
# mu1<- mu1[1:n]
M<- rep(c(1:K),each=n/K)

mu1[1:(n/K),]<-1
# mu1[1:(n/K),2]<-2

mu1[(n/K+1):(n/K*2),]<-2
# mu1[(n/K+1):(n/K*2),2]<-1

mu1[(n/K*2+1):(n/K*3),]<-3
# mu1[(n/K*2+1):(n/K*3),2]<-3

# mu1[1:50,]<- 1.5

###

Zmu<- matrix(mu1,n,d)
X<- matrix( rnorm(n*d, c(Zmu), sqrt(tau)),n,d)
V<- matrix( rnorm(p*d),d,p)

Y<- matrix(rnorm(n*p, c(X%*%V), sd= sqrt(sigma2)), n, p)

# plot(mu1)
plot(X,col=M)

###############################




mclust<-Mclust(Y,G = K,modelNames = "VEI")
raw_M<-apply(mclust$z, 1, function(x){c(1:K)[x==max(x)]})


pca_M<-PCAC(Y,K,d=2)
arc_M<- rDARC(Y,d = 2,k = K,meansIni = matrix(1,2),steps = 1E3,sigma2_ini = 0.1,randomStart = F,ver = 2)
# arc_M<- ARC(Y,K,d = 2)

pdf("simX.pdf",6,6)
plot(X,col=M,xlab="X1",ylab = "X2")
dev.off()

pdf("simPCAX.pdf",6,6)
plot(pca_M$X,col=pca_M$M+1,xlab="X1",ylab = "X2")
dev.off()

pdf("simARCX.pdf",6,6)
plot(arc_M$X,col=arc_M$M+1,xlab="X1",ylab = "X2")
dev.off()

# plot(arc_M$mu,type = "p",lwd=10)

adjustedRandIndex(raw_M,M)
adjustedRandIndex(pca_M$M, M)
adjustedRandIndex(arc_M$M, M)


