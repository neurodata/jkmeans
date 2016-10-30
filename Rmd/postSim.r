require("ggplot2")

setwd("~/git/jkmeans/Rmd/")

#good
load("sims/raw_result_2_0_0.5_1_0.Rda")

nseries<- length(raw_test_results)
K<- 2

#likelihood

loglik_GMM<- function(y, mu,sigma2,w){
  diff<- outer(y,mu,"-")
  sum(log(exp(-diff^2/sigma2/2)/sqrt(sigma2)) %*% (as.matrix(w)))
}


#get Bias Variance for J=j in GMM
getBVMean <- function(j){
  #histogram of the mean
  
  muHat<- raw_n$GMMList[[j]]$mu
  wHat<- raw_n$GMMList[[j]]$w
  
  sigma2Hat<- raw_n$GMMList[[j]]$sigma2
  
  cutoff<- 0.05
  exclude<- (wHat[1,]>(1-cutoff) | wHat[1,]<cutoff)
  
  muHat<- muHat[,,!exclude]
  wHat<- wHat[,!exclude]
  sigma2Hat<- sigma2Hat[!exclude]
  yBatch<- raw_n$yBatch[,,!exclude]
  
  batchN<- batchN- sum(exclude)
  
  loglik<- sapply(c(1:batchN), function(i){
    y<- yBatch[,i]
    mu<- muHat[,i]
    w<- wHat[,i]
    sigma2<- sigma2Hat[i]
    loglik_GMM(y,mu,sigma2,w)
    })
  
  meanLoglik<- mean(loglik)
  
  # dfPlot <- data.frame(
  #   "Mu" = c(mu), 
  #   "Component" = as.factor(rep(c(1:K), batchN))
  # )
  # 
  # meanHist<- ggplot(dfPlot, aes(x=Mu, fill=Component)) + geom_histogram(alpha=0.5, position="identity", bins = 100) + theme_bw()
  
  #bias variance of the mean
  
  Emuhat<- rowMeans(muHat)
  bias <- Emuhat - raw_n$mu0
  variance <-  rowMeans((muHat - Emuhat)^2)
  
  MSE <- rowMeans((muHat - raw_n$mu0)^2)
  
  
  # bias variance of the zeta
  # zeta<- raw_n$GMMList[[J]]$zeta
  # 
  # #weighted Y hat
  # Yhat<- t(sapply(c(1:n) ,function(i){
  #    (t(zeta[i,,]) %*% mu[,,i])
  # }))
  # 
  # EYhat<- rowMeans(Yhat)
  # trueYMean <- sliceMeans(yBatch)
  # variance<-  rowMeans((Yhat - EYhat)^2)
  # bias <-  EYhat- trueYMean
  # # BV<-   rowMeans((Yhat - c(trueYMean))^2)
  # # bias^2 + variance -BV
  
  list("variance"=variance, "bias"=bias, "MSE"=MSE, "meanLoglik"= meanLoglik)
}


n<-1
raw_n<- raw_test_results[[n]]

hist(raw_n$M)

batchN<- dim(raw_n$yBatch)[3]
yBatch <- raw_n$yBatch

sliceMeans <- function(x){
  matrix( colMeans(matrix(c(x),nrow = dim(x)[3],byrow = T)), nrow = dim(x)[1])
}


getBVMean(1)
getBVMean(2)







