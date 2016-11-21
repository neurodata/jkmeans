#

require("sparcl")

n<- 300
p<- 2000
d<- 3
sigma2<- 2
tau<- 0.5
Sigma<- diag(tau,d)

###################################

K<- 3
mu1<- rep(K,n)
# randomMu<- unlist(sapply(1:K, function(x){rep(x, runif(1,n/10,n))}))
# mu1[1:length(randomMu)]<- randomMu
# mu1<- mu1[1:n]
mu1[1:(n/K)]<-1
mu1[(n/K+1):(n/K*2)]<-2
mu1[(n/K*2+1):(n/K*3)]<-3

###############################


Zmu<- matrix(mu1,n,d)
X<- matrix( rnorm(n*d, c(Zmu), sqrt(tau)),n,d)
V<- matrix( rnorm(p*d),d,p)

Y<- matrix(rnorm(n*p, c(X%*%V), sd= sqrt(sigma2)), n, p)

# R<- svd(matrix(rnorm(n^2),n))$u

# Y<- R%*%Y
X0<- X

###############
#Simple PCA+ Cluster



#######

####initialization####
setwd("~/git/jkmeans/test/")
source("kmpp.r")
require("jkmeans")
require("mclust")
source("ARC.R")
source("PCAC.R")

pca_M<-PCAC(Y,K,d)
arc_M<-ARC(Y,K,d,randomStart = T)
pcaGMM_M<-PCAC(Y,K,d,repul = F)

ldgmm_M<-ARC(Y,K,d,randomStart = T,repul = F)


pdf("originalY2d.pdf",6,6)
plot(Y[,1:2],col=mu1, xlab="",ylab="")
dev.off()

pdf("originalX.pdf",6,6)
plot(X,col=mu1, xlab="",ylab="")
dev.off()
# plot(Y[,1:2],col=mu1)

pdf("pca_gmm.pdf",6,6)
plot(pcaGMM_M$X, col=pcaGMM_M$M +1, xlab="",ylab="")
dev.off()

pdf("pca_arc.pdf",6,6)
plot(pca_M$X, col=pca_M$M +1, xlab="",ylab="")
dev.off()

pdf("ld_gmm.pdf",6,6)
plot(ldgmm_M$X, col=ldgmm_M$M +1, xlab="",ylab="")
dev.off()

pdf("ld_arc.pdf",6,6)
plot(arc_M$X, col=arc_M$M +1, xlab="",ylab="")
dev.off()


pdf("ld_arc_truth.pdf",6,6)
plot(arc_M$X, col=4-mu1, xlab="",ylab="")
dev.off()


adjustedRandIndex(pca_M$M,mu1)
adjustedRandIndex(pcaGMM_M$M,mu1)
adjustedRandIndex(ldgmm_M$M,mu1)
adjustedRandIndex(arc_M$M,mu1)


plot(pca_M$M)
plot(arc_M$M)


data("iris")

Y<- as.matrix(iris[,1:4])
# iris[,5]

# plot(X)
K<-3
ini<- as.matrix(kmpp(Y,K))
raw_Clust<- jkmeansEM(Y,k = K,j= 2,1000,tol = 1E-15,fixW = F, meansIni = ini ,useKmeansIni = F,sigma2_ini = 0.1,normalizeZeta = T)
raw_M<- raw_Clust$M

mclust<-Mclust(Y,G = K)
raw_M<-apply(mclust$z, 1, function(x){c(1:K)[x==max(x)]})
####

pca_M<-PCAC(Y,K)
arc_M<-ARC(Y,K)
###


adjustedRandIndex(raw_M, iris[,5])
adjustedRandIndex(pca_M$M, iris[,5])
adjustedRandIndex(arc_M$M, iris[,5])

# plot()
sum((3-arc_M$M)!=as.numeric(iris[,5]))
sum((3-pca_M$M)!=as.numeric(iris[,5]))
sum((raw_M)!=as.numeric(iris[,5]))

# plot(arc_M)
# plot(pca_M)
# plot(raw_M)


#sparse Cl
km.perm <- KMeansSparseCluster.permute(Y,K,wbounds=seq(3,7,len=15),nperms=5)
km.out <- KMeansSparseCluster(Y,K,wbounds=km.perm$bestw)
sparse_M<- as.numeric(km.out[[1]]$Cs)
  
adjustedRandIndex(sparse_M, iris[,5])

plot(sparse_M)
