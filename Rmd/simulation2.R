simulation_raw<- function(K=2, balanced = T, sigma = 1, initial_truth = T, fix_sigma2 = F){
  
  experiment<- function(n){
    
    p<- 1
    batchN<- 1000
    
    yBatch<- array(0,dim = c(n*K,p,batchN))
    
    M <- numeric() #label
    
    for(k in 1:K){
      if(balanced)
        M<- c(M, rep(k,n))
      else{
        if(k==K & K%%2==1){
          M<- c(M, rep(k, n))
        }else{
          M<- c(M, rep(k, n/2 * ((k %% 2)*2+1)  ))
        }
      }
    }
    
    w<- table(M)/sum(table(M))
    
    mu<- 1*M
    mu0<- 1*c(1:K)
    
    for(b in 1:batchN){
      yBatch[,,b]<- matrix( rnorm(n*p*K,mu,sigma), n*K, p)
    }
    
    useKmeansIni =  !initial_truth
    
    kMeansList <- lapply(c(1:K), function(J){
      jkmeans::jkmeansEMBatch(y= yBatch, k=K, j = J, steps =  1000,tol = 1E-10,useKmeansIni = useKmeansIni, meansIni = as.matrix(mu0), fixW = T,sigma2_ini = sigma^2, fixSigma2 = fix_sigma2,flexJ = F,zetaTrunc = 0.1)
    })
    
    GMMList <- lapply(c(1:K), function(J){
      jkmeans::jkmeansEMBatch(y = yBatch, k=K, j = J, steps =  1000,tol = 1E-10,useKmeansIni = useKmeansIni, meansIni = as.matrix(mu0), fixW = F,sigma2_ini = sigma^2, fixSigma2 = fix_sigma2,flexJ = F,zetaTrunc = 0.1)
    })
    
    list("M"= M,"mu0"=mu0, "yBatch"= yBatch, "kMeansList"=kMeansList, "GMMList"= GMMList)
    
  }
  
  nSeries<- c(seq(20,200,by = 20),300,400,500)
  
  raw_test_results <- lapply(nSeries, experiment)
  
  filename = paste("sims/raw_result_",paste(K,balanced ,sigma , initial_truth, as.numeric(fix_sigma2),sep =  "_"),".Rda",sep="")
  save(raw_test_results, file=filename)
  
  # result
}
