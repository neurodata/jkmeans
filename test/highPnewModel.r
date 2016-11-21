require("jkmeans")

n<- 100
K<- 2
p<- 100

sigma<- 1
mu1<- rep(K,n*K)


randomMu<- unlist(sapply(1:K, function(x){rep(x, runif(1,n/10,n))}))
mu1[1:length(randomMu)]<- randomMu


kmpp <- function(X, k) {
  n <- nrow(X)
  C <- numeric(k)
  C[1] <- sample(1:n, 1)
  
  for (i in 2:k) {
    dm <- abs(outer(X, X[C, ], "-"))
    pr <- apply(dm, 1, min)
    pr[C] <- 0
    C[i] <- sample(1:n, 1, prob = pr)
  }
  
  X[C, ]
}

mu<- matrix(0,n*K,p)
mu[,1:5]<- mu1

mu[10:20,1:5]<- 1.5

y<- matrix( rnorm(n*K*p, mu,sd = (sigma)), n*K)
# hist(y,breaks = 100)


svdY<- svd(y)

plot(svdY$d)

V<- svdY$v[1,]

X<- as.matrix(svdY$u[,1]*svdY$d[1])

k<- 3

ini<- as.matrix(kmpp(X,k))

trace_mu<- numeric()

for(i in 1:100){
  jk<- jkmeansEM(X,k = k,j= 1,1000,tol = 1E-15,fixW = F, meansIni = ini ,useKmeansIni = F,sigma2_ini = 0.1,normalizeZeta = T)
  
  Zmu<- jk$zeta %*% jk$mu
  V<- solve(t(Zmu)%*%Zmu ,t(Zmu)%*%y)
  X<- t(solve(V%*%t(V), t(y%*%t(V))))
  ini <- as.matrix(jk$mu[order(jk$mu[,1]),])
  
  trace_mu<- rbind(trace_mu, t(ini))
  print(min(sum(jk$M+1 != mu1)/n/K, 1-sum(jk$M+1 != mu1)/n/K))
}

est<- jk$mu%*%(V)
est[,1]

jk$w

ts.plot(trace_mu)

# require("mclust")

# adjustedRandIndex(jk$M,mu1)

