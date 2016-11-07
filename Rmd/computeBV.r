require("ggplot2")

setwd("~/git/jkmeans/Rmd/")

#good
load("sims/raw_result_2_0_0.5_1_0.Rda")

# load("sims/raw_result_2_0_0.5_0_1.Rda")

nseries<- length(raw_test_results)
K<- 2

#likelihood

loglik_GMM<- function(y, mu,sigma2,w){
  diff<- outer(y,mu,"-")
  sum(log(exp(-diff^2/sigma2/2)/sqrt(sigma2)) %*% (as.matrix(w)))
}


#get Bias Variance for J=j in GMM
getBVMean <- function(j,n){
  #histogram of the mean
  
  raw_n<- raw_test_results[[n]]
  
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
  
  # loglik<- sapply(c(1:batchN), function(i){
  #   y<- yBatch[,i]
  #   mu<- muHat[,i]
  #   w<- wHat[,i]
  #   sigma2<- sigma2Hat[i]
  #   loglik_GMM(y,mu,sigma2,w)
  #   })
  # 
  # meanLoglik<- mean(loglik)
  # 
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
  
  data.frame("variance"=variance, "bias"=bias, "MSE"=MSE, "sigma2"= mean(sigma2Hat)) #, "meanLoglik"= meanLoglik)
}

batchN<- 1000


res<- lapply(c(1:K), function(j){
  lapply(c(1:nseries), function(n){
    getBVMean(j,n)
  })
})

nSeries<- c(seq(20,200,by = 20),300,400,500)


df<- data.frame()

for(j in 1:K){
  variance<- sapply(res[[j]],function(x){x$variance})
  bias<- sapply(res[[j]],function(x){x$bias})
  mse<- sapply(res[[j]],function(x){x$MSE})
  sigma2<- sapply(res[[j]],function(x){x$sigma2})
  
  df<- rbind( df, data.frame("Quantity"="Variance", "Value"=variance[1,] ,"n"=nSeries, "J"= as.factor(j)))
  
  df<- rbind( df, data.frame("Quantity"="Bias2", "Value"=abs(bias[1,])^2 ,"n"=nSeries, "J"= as.factor(j)))
  df<- rbind( df, data.frame("Quantity"="MSE", "Value"=mse[1,] ,"n"=nSeries, "J"= as.factor(j)))
  df<- rbind( df, data.frame("Quantity"="Sigma2", "Value"=sigma2[1,] ,"n"=nSeries, "J"= as.factor(j)))
  
}

pdf("BV2wSigma2.pdf")
p<- ggplot(data=df, aes())
p<- p+  geom_line(aes(x= n, y=Value, colour=J,linetype=Quantity )) +facet_wrap(~J)
p+theme_bw()
dev.off()


pdf("BV2.pdf",6,3)
p<- ggplot(data=df[df$Quantity!="Sigma2",], aes())
p<- p+  geom_line(aes(x= n, y=Value, colour=J,linetype=Quantity )) +facet_wrap(~J)
p+theme_bw()
dev.off()

pdf("BV2Overlay.pdf")
p<- ggplot(data=df[df$Quantity!="Sigma2",], aes())
p<- p+  geom_line(aes(x= n, y=Value, colour=J,linetype=Quantity ))
p+theme_bw()
dev.off()




raw_n<- raw_test_results[[1]]
muHat<- raw_n$GMMList[[1]]$mu
mus1 <- t(sapply(c(1:batchN), function(i){muHat[,,i]}))


muHat<- raw_n$GMMList[[2]]$mu
mus2 <- t(sapply(c(1:batchN), function(i){muHat[,,i]}))

df1<- data.frame(melt(mus1),"J"=1)
df2<- data.frame(melt(mus2),"J"=2)

df<- rbind(df1,df2)

names(df)[2]<- "mu"
df$mu<- as.factor(df$mu)
df$J<- as.factor(df$J)

pdf("BVMeanDist.pdf",6,3)
p<- ggplot(data=df, aes())
p+  geom_histogram(aes(x=value,fill=mu),bins=100, alpha = 0.8,position = 'identity')+facet_wrap(~J)
dev.off()



raw_n<- raw_test_results[[2]]
muHat<- raw_n$GMMList[[1]]$mu
mus1 <- t(sapply(c(1:batchN), function(i){muHat[,,i]}))


muHat<- raw_n$GMMList[[2]]$mu
mus2 <- t(sapply(c(1:batchN), function(i){muHat[,,i]}))

df1<- data.frame(melt(mus1),"J"=1)
df2<- data.frame(melt(mus2),"J"=2)

df<- rbind(df1,df2)

names(df)[2]<- "mu"
df$mu<- as.factor(df$mu)
df$J<- as.factor(df$J)

pdf("BVMeanDist2.pdf",6,3)
p<- ggplot(data=df, aes())
p+  geom_histogram(aes(x=value,fill=mu),bins=100, alpha = 0.8)+facet_wrap(~J)
dev.off()



raw_n<- raw_test_results[[3]]
muHat<- raw_n$GMMList[[1]]$mu
mus1 <- t(sapply(c(1:batchN), function(i){muHat[,,i]}))


muHat<- raw_n$GMMList[[2]]$mu
mus2 <- t(sapply(c(1:batchN), function(i){muHat[,,i]}))

df1<- data.frame(melt(mus1),"J"=1)
df2<- data.frame(melt(mus2),"J"=2)

df<- rbind(df1,df2)

names(df)[2]<- "mu"
df$mu<- as.factor(df$mu)
df$J<- as.factor(df$J)

pdf("BVMeanDist3.pdf",6,3)
p<- ggplot(data=df, aes())
p+  geom_histogram(aes(x=value,fill=mu),bins=100, alpha = 0.8)+facet_wrap(~J)
dev.off()



raw_n<- raw_test_results[[6]]
muHat<- raw_n$GMMList[[1]]$mu
mus1 <- t(sapply(c(1:batchN), function(i){muHat[,,i]}))


muHat<- raw_n$GMMList[[2]]$mu
mus2 <- t(sapply(c(1:batchN), function(i){muHat[,,i]}))

df1<- data.frame(melt(mus1),"J"=1)
df2<- data.frame(melt(mus2),"J"=2)

df<- rbind(df1,df2)

names(df)[2]<- "mu"
df$mu<- as.factor(df$mu)
df$J<- as.factor(df$J)

pdf("BVMeanDist6.pdf",6,3)
p<- ggplot(data=df, aes())
p+  geom_histogram(aes(x=value,fill=mu),bins=100, alpha = 0.8)+facet_wrap(~J)
dev.off()

