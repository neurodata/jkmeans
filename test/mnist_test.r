# .rs.restartR()
setwd("~/git/jkmeans/test/")
source("kmpp.r")
require("jkmeans")
require("mclust")



mnist<- read.csv("../mnist_data/mnist.csv",header = F)

# pick<- c(0,3,7,9)
pick<- c(0:3)
p<- 28^2


# subset<- c(sapply(pick*1000, function(x){x+ c(1:1000)}))

n<-100
s<- sample(1:1000,n,replace = F)
# s<- c(1:1000)

subset<- c(sapply(pick*1000, function(x){x+ s}))
mnist_subset<- mnist[,subset]
Y<-t(mnist_subset)
Y<- Y/255
K<- length(pick)

jk<- jkmeansEM(Y,k = K,useKmeansIni = T,meansIni = matrix(rep(1,K*p),nrow = K),sigma2_ini = 0.1,fixW = F,tau = 1000)

gmm<- jkmeansEM(Y,k = K,useKmeansIni = T,meansIni = matrix(rnorm(K*p),nrow = K),sigma2_ini = 0.1,fixW = F,tau = 1,steps = 1000)

kmeans <- jkmeansEM(Y,k = K,useKmeansIni = T,meansIni = matrix(rep(1,K*p),nrow = K),sigma2_ini = 0.1,fixW = T,tau = 1000)


M<- rep(pick,each=n)

adjustedRandIndex(jk$M,M)
adjustedRandIndex(gmm$M,M)
adjustedRandIndex(kmeans$M,M)

plot(jk$M)

# # mclust<-Mclust(Y,G = K)
# 
# 
# testMNIST<- function(n=100){
#   
# 
#   Y= Y+rnorm(length(Y),mean = 0,sd = 20)
#   Y= Y/255
#   M<- rep(pick,each=n)
#   K<- length(pick)
#   
#   # Y<- (Y-mean(Y))/sd(Y)
#   
#   adjustedRandIndex(raw$M,M)
#   print("raw") 
#   pca<- PCAC(Y,K,d = 6)
#   adjustedRandIndex(pca$M, M)
#   print("pca") 
#   arcM<- rDARC(Y,d = 6,k = K,meansIni = matrix(1,2),steps = 1E3,sigma2_ini = 0.1,randomStart = F,ver = 1,useEstep = F)
#   arc<- rDARC(Y,d = 6,k = K,meansIni = matrix(1,2),steps = 1E3,sigma2_ini = 0.1,randomStart = F,ver = 1,useEstep = T)
#   
#   adjustedRandIndex(arc$M,M)
#   
#   print("arc") 
#   
#   c(adjustedRandIndex(raw$M,M),
#     adjustedRandIndex(pca$M, M),
#     adjustedRandIndex(arcM$M, M),
#     adjustedRandIndex(arc$M, M))
# }
# 
# t50<- sapply(c(1:10), function(x)testMNIST(100))
# 
# 
# drowSDs<-function(x){
#   apply(x, 1, function(x)sd(x,na.rm = T))
# }
# 
# 
# save(t50, file="t50.Rda")
# 
# load(file="t50.Rda")
# 
# t50[t50==0]<- NA
# rowMeans(t50,na.rm = T)
# rowSDs(t50)
# 
# 
# 
# 
# 
# require("sparcl")
# 
# #sparse Cl
# km.perm <- KMeansSparseCluster.permute(Y,K,wbounds=seq(3,7,len=15),nperms=5)
# km.out <- KMeansSparseCluster(Y,K,wbounds=km.perm$bestw)
# sparse_M<- as.numeric(km.out[[1]]$Cs)
# 
# adjustedRandIndex(sparse_M, M)
