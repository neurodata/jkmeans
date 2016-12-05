require('ggplot2')

sigma<- 0.5
mu<- c(rep(1,100),rep(2,100))
Y<- rnorm(200,mu,sigma)


# dataHist<- data.frame("Center"=as.factor(mu),    "y"= Y)
# ggplot(dataHist, aes(y, fill = Center)) + geom_histogram(alpha = 0.5,bins = 30,  position="identity")

K<-2



sigma<- 0.3
repMu3<-sapply(c(1:2000), function(x){
  Y<- rnorm(200,mu,sigma)
  raw<- jkmeansEM(as.matrix(Y),k = K,useKmeansIni = T,meansIni = matrix(rep(1,2),nrow = K),sigma2_ini = 0.1,fixW = F)
  raw$mu
})

gmmMu3<-sapply(c(1:2000), function(x){
  Y<- rnorm(200,mu,sigma)
  fi<- Mclust(as.matrix(Y),2,modelNames="E")
  fi$parameters$mean
})



sigma<- 0.5
repMu5<-sapply(c(1:2000), function(x){
  Y<- rnorm(200,mu,sigma)
  raw<- jkmeansEM(as.matrix(Y),k = K,useKmeansIni = T,meansIni = matrix(rep(1,2),nrow = K),sigma2_ini = 0.1,fixW = F)
  raw$mu
})


gmmMu5<-sapply(c(1:2000), function(x){
  Y<- rnorm(200,mu,sigma)
  fi<- Mclust(Y,2,modelNames="E")
  fi$parameters$mean
})

sigma<- 1
repMu10<-sapply(c(1:2000), function(x){
  Y<- rnorm(150,mu,sigma)
  raw<- jkmeansEM(as.matrix(Y),k = K,useKmeansIni = T,meansIni = matrix(rep(1,2),nrow = K),sigma2_ini = 0.1,fixW = F)
  raw$mu
})


gmmMu10<-sapply(c(1:2000), function(x){
  Y<- rnorm(200,mu,sigma)
  fi<- Mclust(Y,2,modelNames="E")
  fi$parameters$mean
})


sigma<- 1.5
repMu15<-sapply(c(1:2000), function(x){
  Y<- rnorm(150,mu,sigma)
  raw<- jkmeansEM(as.matrix(Y),k = K,useKmeansIni = T,meansIni = matrix(rep(1,2),nrow = K),sigma2_ini = 0.1,fixW = F)
  raw$mu
})


gmmMu15<-sapply(c(1:2000), function(x){
  Y<- rnorm(200,mu,sigma)
  fi<- Mclust(Y,2,modelNames="E")
  fi$parameters$mean
})


dataHist<- data.frame("Component"=as.factor(rep(c(1,2),2000*3)),  
                      "Data"= c(c(rnorm(2000,mu,0.3)),c(rnorm(2000,mu,0.5)),c(rnorm(2000,mu,1))), "Std.Dev"=as.factor(rep(c(0.3,0.5,1),each=4000)))
yPlot<- ggplot(dataHist, aes(Data, fill = Component)) + geom_histogram(alpha = 0.5,bins = 40,  position="identity")+xlim(c(-1,4))+ geom_vline(xintercept = 1,lty=2)+ geom_vline(xintercept = 2,lty=2)+theme_bw() + facet_grid(. ~ Std.Dev)



dataHist<- data.frame("Component"=as.factor(rep(c(1,2),2000*3)),   "Mean Estimate"= c(c(repMu3),c(repMu5),c(repMu10)), "Std.Dev"=as.factor(rep(c(0.3,0.5,1),each=4000)))
mPlot<- ggplot(dataHist, aes(Mean.Estimate, fill = Component)) + geom_histogram(alpha = 0.5,bins = 50,  position="identity")+xlim(c(-1,4))+ geom_vline(xintercept = 1,lty=2)+ geom_vline(xintercept = 2,lty=2)+theme_bw() + facet_grid(. ~ Std.Dev)



dataHist<- data.frame("Component"=as.factor(rep(c(1,2),2000*3)),   "Mean Estimate"= c(c(gmmMu3),c(gmmMu5),c(gmmMu10)), "Std.Dev"=as.factor(rep(c(0.3,0.5,1),each=4000)))
gmmPlot<- ggplot(dataHist, aes(Mean.Estimate, fill = Component)) + geom_histogram(alpha = 0.5,bins = 50,  position="identity")+xlim(c(-1,4))+ geom_vline(xintercept = 1,lty=2)+ geom_vline(xintercept = 2,lty=2)+theme_bw() + facet_grid(. ~ Std.Dev)


pdf("repulData.pdf",12,4)
yPlot
dev.off()

pdf("repulMean.pdf",12,4)
mPlot
dev.off()


pdf("gmmMean.pdf",12,4)
gmmPlot
dev.off()

mu<- rep(c(1,2),1000)


