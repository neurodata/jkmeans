    ###Prerequisite
    ###Install the package if you haven't
    setwd("~/git/")
    require('devtools')
    build('jkmeans')
    install.packages("jkmeans_1.0.tar.gz", repos = NULL, type = "source")

### 1. Generate Data

Let's generate data:

``` r
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
```

Here is a view of the two clusters

``` r
dataHist<- data.frame("Mu"=as.factor(mu),    "y"= yBatch[,,1])

ggplot(dataHist, aes(y, fill = as.factor(Mu))) + geom_histogram(alpha = 0.2,bins = 30,  position="identity")
```

![](convergenceRatesBigK3_files/figure-markdown_github/unnamed-chunk-3-1.png)

functions to compute misclassification error and rmse of mean estimate: NB: GMM can get stuck in saddle point at (w1=0, w2=1), so I removed the test result when it happened.

``` r
computeMError<- function(jk, n,batchN){
  error<- numeric(batchN)
  for(b in 1:batchN){
    tM0<- rank(jk$mu[,1,b])-1
    trueM<- rep(tM0,each=n)
    error[b]<- sum ( jk$M[,b] != trueM) /K/n
    if( min(jk$w[,b])<0.01){
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
    if( min(jk$w[,b])<0.01)
      rmse[b]<- NA
  }
  rmse
}
```

### 2. Test jk-means/ j-sparse-GMM with different n's

``` r
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
```

``` r
#code to run on an increasing series of n
#ran on cluster

nSeries<- c(seq(20,500,by = 20))

if(FALSE)
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
  
  save(result, file="resultRatesBigK3.Rda")
}
```

Misclassification error
=======================

plot without & with pointwise 95% confidence band

``` r
load("resultRatesBigK3.Rda")

plot1<- data.frame( "n" = rep(nSeries, 4),
            "MC14Error"= c(result$MCEmean),
            "MC14ErrorL"= c(result$MCEq25),
            "MC14ErrorU"= c(result$MCEq975),
            "Model"= rep(c("K-Means","W-constraint-GMM (22-Means)","1-sparse-GMM","GMM"),each= length(nSeries))
            )

p<- ggplot(data=plot1, aes())
p<- p+  geom_line(aes(x= n, y=MC14Error, colour=Model)) 
pE<- p+  geom_line(aes(x= n, y=MC14ErrorL, colour=Model),linetype=2) 
pE<- pE+  geom_line(aes(x= n, y=MC14ErrorU, colour=Model),linetype=2) 

# p+ geom_errorbar( aes(x=n,ymax = MC14ErrorU, ymin=MC14ErrorL, colour=Model),width=0.2)
p+ theme_bw()
```

![](convergenceRatesBigK3_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
pE + theme_bw()
```

![](convergenceRatesBigK3_files/figure-markdown_github/unnamed-chunk-7-2.png)

RMSE for the mean estimate
==========================

plot without & with pointwise 95% confidence band

``` r
plot1<- data.frame( "n" = rep(nSeries, 4),
            "RMSE"= c(result$RMSEmean),
            "RMSEL"= c(result$RMSEq25),
            "RMSEU"= c(result$RMSEq975),
            "Model"= rep(c("K-Means","W-constraint-GMM (22-Means)","1-sparse-GMM","GMM"),each= length(nSeries))
            )

p<- ggplot(data=plot1, aes())
p<- p+  geom_line(aes(x= n, y=RMSE, colour=Model)) 

pE<- p+  geom_line(aes(x= n, y=RMSEL, colour=Model),linetype=2) 
pE<- pE+  geom_line(aes(x= n, y=RMSEU, colour=Model),linetype=2) 


p+ theme_bw()
```

![](convergenceRatesBigK3_files/figure-markdown_github/unnamed-chunk-8-1.png)

``` r
pE + theme_bw()
```

![](convergenceRatesBigK3_files/figure-markdown_github/unnamed-chunk-8-2.png)