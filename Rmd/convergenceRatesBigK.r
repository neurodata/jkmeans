
require('ggplot2')
require('reshape')
require('jkmeans')

setwd("~/git/jkmeans/Rmd/")

n<- 100
K<- 10
p<- 1
batchN<- 1000

yBatch<- array(0,dim = c(n*K,p,batchN))

mu0<- c(1:K)
mu<- matrix( 0,n*K,p)
mu <- rep(mu0, each = n)
sigma<-  0.5

for(b in 1:batchN){
  yBatch[,,b]<- matrix( rnorm(n*p*K,mu,sigma), n*K, p)
}

mu_ini<- matrix(mu0,K,p)



dataHist<- data.frame("Mu"=as.factor(mu),    "y"= yBatch[,,1])

ggplot(dataHist, aes(y, fill = as.factor(Mu))) + geom_histogram(alpha = 0.2,bins = 30,  position="identity")


computeMError<- function(jk, n,batchN){
  error<- numeric(batchN)
  for(b in 1:batchN){
    tM0<- rank(jk$mu[,1,b])-1
    trueM<- rep(tM0,each=n)
    error[b]<- sum ( jk$M[,b] != trueM) /K/n
    if( min(jk$w[,b])<0.1){
        error[b]<- NA
    }
  }
  error
}

computeRMSE<- function(jk, mu0, batchN){
  rmse<- numeric(batchN)
  for(b in 1:batchN){
    newOrder<- order(jk$mu[,1,b],decreasing = F)
    
    muEst<-as.matrix(jk$mu[,,b])[newOrder,]
   
    rmse[b]<- sqrt( mean((muEst - mu0)^2))
    if( min(jk$w[,b])<0.1)
      rmse[b]<- NA
  }
  rmse
}

experiment<- function(n){
  
  p<- 1
  batchN<- 1000
  
  yBatch<- array(0,dim = c(n*K,p,batchN))
  K<- 10
  
  mu0<- c(1:K)
  mu<- matrix( 0,n*K,p)
  mu <- rep(mu0, each = n)
  sigma<-  0.5


  for(b in 1:batchN){
    yBatch[,,b]<- matrix( rnorm(n*p*K,mu,sigma), n*K, p)
  }
  
  mu_ini<- matrix(mu0,K,p)
  
  
  jk12 <- jkmeans::jkmeansEMBatch(yBatch, k=K, j = 1,  1000,tol = 1E-10,useKmeansIni = F, meansIni = mu_ini, fixW = T)
  jk22 <- jkmeans::jkmeansEMBatch(yBatch, k=K, j = K,  1000,tol = 1E-10,useKmeansIni = F, meansIni = mu_ini, fixW = T)
  jGMM12 <- jkmeans::jkmeansEMBatch(yBatch, k=K, j = 1,  1000,tol = 1E-10,useKmeansIni = F, meansIni = mu_ini, fixW = F)
  jGMM22 <- jkmeans::jkmeansEMBatch(yBatch, k=K, j = K,  1000,tol = 1E-10,useKmeansIni = F, meansIni = mu_ini, fixW = F)

  
  error <- cbind( computeMError(jk12,n, batchN), computeMError(jk22,n,batchN), computeMError(jGMM12,n,batchN),computeMError(jGMM22,n,batchN))
  rmse <- cbind( computeRMSE(jk12,mu0,batchN), computeRMSE(jk22,mu0,batchN), computeRMSE(jGMM12,mu0,batchN),computeRMSE(jGMM22,mu0,batchN))
    
  list("MCE"=error, "RMSE"=rmse)
  
}
  


nSeries<- c(seq(20,500,by = 20))

{
  
  MCEmean<- numeric()
  MCEq25<- numeric()
  MCEq975<- numeric()
  RMSEmean<- numeric()
  RMSEq25<- numeric()
  RMSEq975<- numeric()
  
  
  for(n in nSeries){
    test<-experiment(n)
    
    MCEmean <- rbind(MCEmean,colMeans(test$MCE,na.rm = T))
    MCEq25 <- rbind(MCEq25, apply(test$MCE, MARGIN = 2, function(x){quantile(x,probs = 0.025,na.rm = T)}))
    MCEq975 <- rbind(MCEq975, apply(test$MCE, MARGIN = 2, function(x){quantile(x,probs = 0.975,na.rm = T)}))
    
    RMSEmean<- rbind(RMSEmean, colMeans(test$RMSE,na.rm = T))
    RMSEq25<- rbind(RMSEq25, apply(test$RMSE, MARGIN = 2, function(x){quantile(x,probs = 0.025,na.rm = T)}))
    RMSEq975<-  rbind(RMSEq975, apply(test$RMSE, MARGIN = 2, function(x){quantile(x,probs = 0.975,na.rm = T)}))
    
    print(n)
  }
  
  
  result<- list("MCEmean"=MCEmean,
                "MCEq25"=MCEq25,
                "MCEq975"=MCEq975,
                "RMSEmean"=RMSEmean,
                "RMSEq25"=RMSEq25,
                "RMSEq975"=RMSEq975)
  
  save(result, file="resultRatesBalanced.Rda")
}