require("pixmap")
require("jkmeans")
require("mclust")

setwd("~/git/jkmeans/lfwcrop_face/lfwcrop_grey/pick/")

faces<- sapply(list.files("."),function(f){
  x=read.pnm(file = f)
  c(x@grey)
})


M<- sapply(list.files("."), function(x){substr(x,1,5)})

X<- t(faces)

nf<- 64

K<- 2

plotFace<- function(x){
  image(matrix(x,64,64),col = grey(seq(0, 1, length = 256)))
}

raw<- jkmeansEM(X,k = K,useKmeansIni = T,meansIni = matrix(rep(1,2*10),nrow = 2),sigma2_ini = 0.1,fixW = F)

plotFace(raw$mu[1,])
plotFace(raw$mu[2,])


plot(raw$M)

adjustedRandIndex(raw$M,M)


arcM<- rDARC(X,d = 9,k = K,meansIni = matrix(1,2),steps = 5E3,sigma2_ini = 0.1,randomStart = F,ver = 1,useEstep = F)

adjustedRandIndex(arcM$M,M)


plotFace(arcM$mu[1,]%*%arcM$EV)
plotFace(arcM$mu[2,]%*%arcM$EV)

plot(arcM$M)



specClust<-function(X1,K,tau=50){
  N<- nrow(X1)
  dist<- matrix(0,N,N)
  for(i in 1:N){
    for(j in 1:i){
      dist[i,j]<- exp(-tau* sum((X1[i,]-X1[j,])^2))
      dist[j,i]<- dist[i,j]
    }
  }
  
  D<- (colSums(dist))
  L<- t(dist/sqrt(D))/sqrt(D)

  svdL<- svd(L,nu = K)
  X2<- svdL$u
  
  nX2<- X2/ sqrt( rowSums(X2^2))
  
  clust<-kmeans(nX2,K)
  clust$cluster
}


spcM<-specClust(X,2, 0.0001)

adjustedRandIndex(spcM,M)

plot(spcM)
plot(as.numeric(as.factor(M)))

miscError<- function(x){
  rate=sum(x != as.numeric(as.factor(M)))/N
  min(rate,1-rate)
}

miscError(raw$M)
miscError(arcM$M)
miscError(spcM)
