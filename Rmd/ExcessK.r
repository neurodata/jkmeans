require("jkmeans")

n<- 100
K<- 5

sigma<- 1
mu<- rep(c(1:K),each=n)
y<- matrix( rnorm(n*K, mu,sd = (sigma)))
hist(y,breaks = 100)

hist(c(sapply(c(1:K),function(x){ rnorm(1000,x, (sigma2/sqrt(n)))})), breaks = 1000)


k<- K

jk<- jkmeansEM(y,k = k,j= 1,1000,tol = 1E-10,fixW = F, meansIni =  matrix(c(c(1:k))),useKmeansIni = F,sigma2_ini = 0.1)

full<- jkmeansEM(y,k = k,j= k,1000,tol = 1E-10,fixW = F, meansIni =  matrix(c(c(1:k))),useKmeansIni = F,sigma2_ini = 0.1)

oracle <- jkmeansEM(y,k = K,j= K,1000,tol = 1E-10,fixW = F, meansIni =  matrix(c(c(1:K))),useKmeansIni = F,sigma2_ini = 0.1)

flex <- jkmeansEM(y,k = k ,j= k,1000,tol = 1E-10,fixW = F, meansIni =  matrix(c(c(1:k))),useKmeansIni = F,sigma2_ini = 0.1,flexJ = T,zetaTrunc = 0.1)

jk<- jkmeansEM(y,k = k,j= 3,1000,tol = 1E-10,fixW = T, meansIni =  matrix(c(c(1:k))),useKmeansIni = F,sigma2_ini = 1)


table(jk$M)
table(full$M)
table(oracle$M)
table(flex$M)


jk$mu[jk$w>0.01]
full$mu[full$w>0.05]
oracle$mu[oracle$w>0.05]
flex$mu[flex$w>0.01]

