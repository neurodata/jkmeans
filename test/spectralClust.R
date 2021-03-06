N<- 1000
X<- matrix(rnorm(N*2),N)

circle<-function(a,b){
 radius<- (X[,1]^2 +X[,2]^2)
 radius>a^2 & radius<b^2
}

X1<- rbind( X[circle(1,1.5),], X[circle(0,0.6),])

plot(X1)

N<- nrow(X1)
dist<- matrix(0,N,N)
for(i in 1:N){
  for(j in 1:i){
    dist[i,j]<- exp(  - 100* sum((X1[i,]-X1[j,])^2))
    dist[j,i]<- dist[i,j]
  }
}

D<- (colSums(dist))
L<- t(dist/sqrt(D))/sqrt(D)
# image(L)

K<- 2

svdL<- svd(L,nu = K)
X2<- svdL$u

nX2<- X2/ sqrt( rowSums(X2^2))

plot(X2)
# plot(nX2)

clust<-kmeans(nX2,K)

plot(X2,col=clust$cluster)


plot(X1,col=clust$cluster)
