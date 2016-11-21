ARC<-function(Y,K,d=2, randomStart= FALSE, repul=TRUE){
  

  
  svdY<-svd(Y,nu = d,nv = d)
  
  if(d>2)
    X<- svdY$u %*%diag(svdY$d[1:d] )
  else
    X<-   svdY$u * svdY$d[1]
  if(randomStart){
  X<- X + matrix(rnorm(length(c(X)),sd = 1),nrow(X))
  }
  
  

  ini<- as.matrix(kmpp(X,K))
  
  # plot(X,X0)
  
  
  jk<- jkmeansEM(X,k = K,j= 1,1000,tol = 1E-15,fixW = F, meansIni = ini ,useKmeansIni = F,sigma2_ini = 0.1,normalizeZeta = T)
  
  Sigma<-  diag(c(jk$Sigma),d)
  
  Zmu<- jk$zeta %*% jk$mu
  ini_M<- jk$M
  ini_X<- X
  
  
  # plot(jk$M)
  
  # plot(X0,col=mu1)
  # plot(X,col=jk$M+1)
  # plot(X0,col=mu1)
  
  
  # e<- sum(jk$M+1!=mu1)/n
  # min(e,1-e)
  
  #################
  for(i in 1:500){
    
    # for(j in 1:10){
    XX<- t(X)%*%X
    XY<- t(X)%*%Y
    XXinv<- solve(XX+diag(100,d))
    
    #compute expectation
    EV<- solve(XX,XY)
    temp<- XXinv%*%XY
    EVV<- p*sigma2*XXinv +  temp %*% t(temp)
    
    #compute sigma2
    sigma2<- (sum(c(Y)*c(Y))-2* sum(c(t(X)) * c(EV%*%t(Y))) + sum(c(t(X)) * c(EVV%*%t(X)))) /n/p
    
    #get max of X
    
    m<- Y%*%t(EV)/sigma2 + Zmu%*%solve(Sigma)
    v<- solve(EVV/sigma2+solve(Sigma))
    X<- m%*%v
    # }
    
    #clustering
    ini<- as.matrix(jk$mu)
    if(repul){j=1}
    else{j=K}
    jk<- jkmeansEM(X,k = K,j= j,1000,tol = 1E-15,fixW = F, meansIni = ini ,useKmeansIni = F,sigma2_ini = 0.1,normalizeZeta = T)
    Zmu<- jk$zeta %*% jk$mu
    
    Sigma<-  diag(c(jk$Sigma),d)
  }
  return(list("M"=jk$M, "X"=X,"V"=EV,"mu"=jk$mu))
}
