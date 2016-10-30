p<- 1
batchN<- 1000
n<- 100
K<- 2
sigma<- 0.5
fix_sigma2 <- FALSE

yBatch<- array(0,dim = c(n*K,p,batchN))
balanced <- FALSE
initial_truth<- TRUE


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


    list("MCE"=cbind(kMeansError,GMMError), "BayesError"=BayesError,"RMSE"=cbind(kMeansrmse,GMMrmse))

  }

  nSeries<- c(seq(20,200,by = 20),300,400,500)

  {

    MCEmean<- numeric()
    MCEq25<- numeric()
    MCEq975<- numeric()

    BayesError<- numeric()


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

      BayesError<- c(BayesError, test$BayesError)
      print(n)
    }


    result<- list("MCEmean"=MCEmean,
                  "MCEq25"=MCEq25,
                  "MCEq975"=MCEq975,
                  "RMSEmean"=RMSEmean,
                  "RMSEq25"=RMSEq25,
                  "RMSEq975"=RMSEq975,
                  "BayesError"= BayesError
    )
    filename = paste("sims/result_",paste(K,balanced ,sigma , initial_truth, as.numeric(fix_sigma2),sep =  "_"),".Rda",sep="")
    save(result, file=filename)
  }
  # result
}
