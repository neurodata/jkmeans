require("pixmap")
require("jkmeans")
require("mclust")

setwd("~/git/jkmeans/lfwcrop_face/lfwcrop_grey/pick/")

faces<- sapply(list.files("."),function(f){
  x=read.pnm(file = f)
  c(x@grey)
})


M<- sapply(list.files("."), function(x){substr(x,1,5)})

X<- t(faces)

nf<- 64

K<- 2

plotFace<- function(x){
  image(matrix(x,64,64),col = grey(seq(0, 1, length = 256)))
}

raw<- jkmeansEM(X,k = K,useKmeansIni = T,meansIni = matrix(rep(1,2*10),nrow = 2),sigma2_ini = 0.1,fixW = F)

plotFace(raw$mu[1,])
plotFace(raw$mu[2,])


plot(raw$M)

adjustedRandIndex(raw$M,M)


arcM<- rDARC(X,d = 9,k = K,meansIni = matrix(1,2),steps = 5E3,sigma2_ini = 0.1,randomStart = T,ver = 1,useEstep = F)

adjustedRandIndex(arcM$M,M)


plotFace(arcM$mu[1,]%*%arcM$EV)
plotFace(arcM$mu[2,]%*%arcM$EV)

plot(arcM$M)
