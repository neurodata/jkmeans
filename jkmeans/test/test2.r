require("jkmeans")

mu0<- c(1,2)
n<- 40
mu<- c(rep(1,n),rep(2,n*2))
y<- rnorm(n*3, mu)

fit<- jkmeansEM(as.matrix(y),k = 2, j=1, meansIni = as.matrix(mu0),steps = 1000,tol = 1E-15, fixW = F, fixSigma2 = F, sigma2_ini = 1,useKmeansIni = T)

fit$mu
fit$w
fit$sigma2


fit2<- jkmeansEM(as.matrix(y),k = 2, j=2, meansIni = as.matrix(mu0),steps = 1000,tol = 1E-15, fixW = F, fixSigma2 = F, sigma2_ini = 0.01,useKmeansIni = T)

fit$mu
fit$w
fit$sigma2

require("MCMCpack")
rdirichlet(1000, rep(1,5))
