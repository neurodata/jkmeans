n<- 2000
p<- 2000
d<- 2
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


mu<- matrix(mu1,n,2)

X<- matrix(rnorm(n*2,mu),n)
X[,2]<- rnorm(n, X[,2],sd=1)

ini<- as.matrix(kmpp(X,K))


jk<- jkmeansEM(X,k = K,j= 1,1000,tol = 1E-15,fixW = F, meansIni = ini ,useKmeansIni = F,sigma2_ini = 0.1,normalizeZeta = T)

jk$Sigma
