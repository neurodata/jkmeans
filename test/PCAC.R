PCAC<-function(Y,K,d=2, repul=TRUE){
  
  svdY<-svd(Y,nu = d,nv = d)
  
  if(d>2)
    X<- svdY$u %*%diag(svdY$d[1:d] )
  else
    X<-   svdY$u * svdY$d[1]
  
  ini<- as.matrix(kmpp(X,K))
  
  # if(repul){j=1}
  # else{j=K}
  # jk<- jkmeansEM(X,k = K,j= j,1000,tol = 1E-15,fixW = F, meansIni = ini ,useKmeansIni = F,sigma2_ini = 0.1,normalizeZeta = T)
  jk<- jkmeansEM(y= X,k = K,steps = 1000,tol = 1E-15,fixW = F, meansIni = ini ,useKmeansIni = F,sigma2_ini = 0.1)
  
  
  return(list("M"=jk$M, "X"=X,"V"=svdY$v,"mu"=jk$mu))
}