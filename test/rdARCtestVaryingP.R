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

LowDTest<-function(p,   p_trunc= 100){
  n<- 300
  d<- 2
  sigma2<- 1
  tau<- 0.2
  Sigma<- diag(tau,d)
  
  ###################################
  
  K<- 3
  mu1<- matrix(K,n,d)
  # randomMu<- unlist(sapply(1:K, function(x){rep(x, runif(1,n/10,n))}))
  # mu1[1:length(randomMu)]<- randomMu
  # mu1<- mu1[1:n]
  M<- rep(c(1:K),each=n/K)
  M<- c( rep(1,n/K/2), rep(2,n/K/2*3), rep(3,n/K))
  
  mu1[1:(n/K/2),]<-1
  mu1[1:(n/K/2),2]<-2
  
  mu1[(n/K/2+1):(n/K*2),]<-2
  mu1[(n/K/2+1):(n/K*2),2]<-1
  
  mu1[(n/K*2+1):(n/K*3),]<-3
  mu1[(n/K*2+1):(n/K*3),2]<-3
  
  Zmu<- matrix(mu1,n,d)
  X<- matrix( rnorm(n*d, c(Zmu), sqrt(tau)),n,d)
  V<- matrix( rnorm(p*d),d,p)
  # 
  Y<- matrix(rnorm(n*p, c(X%*%V), sd= sqrt(sigma2)), n, p)
  
  
  
  # p_trunc<- p
  Y[,p_trunc:p]<- rnorm(length(Y[,p_trunc:p]),mean=0,sd= sqrt(sigma2))
  
  # Zmu1[,p_trunc:p]<- 0
  # Y<- matrix( rnorm(n*p,c(Zmu1),sd= sqrt(sigma2)), n, p)
  
  
  # plot(mu1)
  
  
  
  
  # mclust<-Mclust(Y[,1:(p_trunc-1)],G = K)
  # raw_M<-apply(mclust$z, 1, function(x){c(1:K)[x==max(x)]})
  
  
  raw<- jkmeansEM(Y,k = K,useKmeansIni = T,meansIni = matrix(rep(1,K*p),nrow = K),sigma2_ini = 0.1,fixW = F)
  raw_M<- raw$M
  
  pca_M<-PCAC(Y,K,d=d)
  
  arc_M<- rDARC(Y,d = d,k = K,meansIni = matrix(1,2),steps = 1E4,sigma2_ini = 0.1,randomStart = F,fixW = F,ver = 1)
  
  # testARC<- lapply(c(1:10),function(x) rDARC(Y,d = 2,k = K,meansIni = matrix(1,2),steps = 1E3,sigma2_ini = 0.1,randomStart = F,ver = 1))
  # l_loglik<- (sapply(testARC, function(x)x$loglik))
  # arc_M<- testARC[[c(1:10)[l_loglik== min(l_loglik)]]]
  
  # plot(X,col=M)
  # plot(pca_M$X,col=pca_M$M+1,xlab="X1",ylab = "X2",lwd=2)
  # plot(arc_M$X,col=arc_M$M+1,xlab="X1",ylab = "X2",lwd=2 )
  #,ylim=range(pca_M$X[,2]))
  
  
  c(adjustedRandIndex(raw_M,M),
    adjustedRandIndex(pca_M$M, M),
    adjustedRandIndex(arc_M$M, M))
}

rowSDs<-function(x){
  apply(x, 1, sd)
}


t50<- sapply(c(1:10), function(x)LowDTest(50,20))

rowMeans(t50)
rowSDs(t50)



t100<- sapply(c(1:10), function(x)LowDTest(100,20))

rowMeans(t100)
rowSDs(t100)


t200<- sapply(c(1:10), function(x)LowDTest(200,20))

rowMeans(t200)
rowSDs(t200)


t500<- sapply(c(1:10), function(x)LowDTest(500,20))

rowMeans(t500)
rowSDs(t500)

t1000<- sapply(c(1:10), function(x)LowDTest(1000,20))

rowMeans(t1000)
rowSDs(t1000)



plot(raw_M)
plot(arc_M$M)


# arc_M$loglik

pdf("simX.pdf",6,6)
plot(X,col=M,xlab="X1",ylab = "X2",lwd=2)
dev.off()
# 
pdf("simPCAX.pdf",6,6)
plot(pca_M$X,col=pca_M$M+1,xlab="X1",ylab = "X2",lwd=2)
dev.off()
# 
pdf("simARCX.pdf",6,6)
plot(arc_M$X,col=arc_M$M+1,xlab="X1",ylab = "X2",lwd=2 )
dev.off()

# plot(arc_M$mu,type = "p",lwd=10)

# arc_M$M

