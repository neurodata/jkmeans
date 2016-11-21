n<- 50
p<- 2000
d<- 2
sigma2<- 3
tau<- 0.05
Sigma<- diag(tau,d)

Zmu<- matrix( rnorm(n*d),n,d)
X<- matrix( rnorm(n*d, c(Zmu), sqrt(tau)),n,d)
V<- matrix( rnorm(p*d),d,p)

Y<- matrix(rnorm(n*p, c(X%*%V), sd= sqrt(sigma2)), n, p)

X0<- X

####

svdY<-svd(Y,nu = d,nv = d)

X<- svdY$u %*%diag(svdY$d[1:d] )
plot(X,X0)
plot(Zmu,X0,xlim=c(-3,3),ylim = c(-3,3))

#################
for(i in 1:1000){
  XX<- t(X)%*%X
  XY<- t(X)%*%Y
  XXinv<- solve(XX+diag(1,d))
  
  #compute expectation
  EV<- solve(XX,XY)
  temp<- XXinv%*%XY
  EVV<- p*sigma2*XXinv +  temp %*% t(temp)
  
  #compute sigma2
  sigma2<- (sum(c(Y)*c(Y))-2* sum(c(t(X)) * c(EV%*%t(Y))) + sum(c(t(X)) * c(EVV%*%t(X)))) /n/p
  
  #get max of X
  
  m<- Y%*%t(EV)/sigma2 + Zmu%*%solve(Sigma)
  v<- solve(EVV/sigma2+solve(Sigma))
  X<- m%*%v
}
###

plot(X,X0,xlim=c(-3,3),ylim = c(-3,3))
