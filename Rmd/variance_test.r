K=2
balanced = F
sigma = 1
initial_truth = T
n= 20

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
    jkmeans::jkmeansEMBatch(yBatch, k=K, j = J,  1000,tol = 1E-10,useKmeansIni = useKmeansIni, meansIni = as.matrix(mu0), fixW = T,sigma2_ini = sigma^2, fixSigma2 = T)
})


GMMList <- lapply(c(1:K), function(J){
  jkmeans::jkmeansEMBatch(yBatch, k=K, j = J,  1000,tol = 1E-10,useKmeansIni = useKmeansIni, meansIni = as.matrix(mu0), fixW = F,sigma2_ini = sigma^2, fixSigma2 = T)
})

tensorMean<- function(x){
  rowMeans( matrix(c(x),nrow = 2))
}

computeVarBias2<- function(x){
  EmuHat<- tensorMean(x)
  variance <- tensorMean((x - EmuHat)^2)
  bias2 <- tensorMean((x - mu0))^2
  return(list("var"=variance, "bias2"=bias2))
}

a<- computeVarBias2(GMMList[[1]]$mu)
b<- computeVarBias2(GMMList[[2]]$mu)

a$bias2
b$bias2

hist(GMMList[[1]]$mu,breaks = 100)
hist(GMMList[[2]]$mu,breaks = 100)

a<- computeVarBias2(kMeansList[[1]]$mu)
b<- computeVarBias2(kMeansList[[2]]$mu)


hist(kMeansList[[1]]$mu,breaks = 100)
hist(kMeansList[[2]]$mu,breaks = 100)

kMeansError<- sapply(kMeansList,function(x){computeMError(x,batchN, M)})
kMeansError<- sapply(kMeansList,function(x){computeMError(x,batchN, M)})

colMeans(kMeansError)

    #functions to compute the error
    computeMError<- function(jk, batchN, trueM){
      error<- numeric(batchN)
      for(b in 1:batchN){
        error[b]<- mean ( jk$M[,b] != (trueM -1 ))
        if( min(jk$w[,b])<0.01){
          error[b]<- NA
        }
      }
      error
    }

    computeRMSE<- function(jk, batchN, trueMu){
      rmse<- numeric(batchN)
      for(b in 1:batchN){
        newOrder<- order(jk$mu[,1,b],decreasing = F)
        muEst<-as.matrix(jk$mu[,,b])[newOrder,]

        rmse[b]<- sqrt( mean((muEst - trueMu)^2))
        if( min(jk$w[,b])<0.01)
          rmse[b]<- NA
      }
      rmse
    }

    #Bayes Error (ugly code, be careful if you want to change)
    diff<- outer( c(yBatch), mu0, "-")
    GaussianLoglik<- -diff^2/sigma^2/2
    loglik<- GaussianLoglik + rep(log(w), each=nrow(diff))
    BayesM<- t(apply(loglik, 1, function(x){c(1:K)[x==max(x)]}))
    BayesError <- mean(BayesM != rep(M,batchN))


    kMeansError<- sapply(kMeansList,function(x){computeMError(x,batchN, M)})
    kMeansrmse<- sapply(kMeansList,function(x){computeRMSE(x,batchN, mu0)})


    GMMError<- sapply(GMMList,function(x){computeMError(x,batchN, M)})
    GMMrmse<- sapply(GMMList,function(x){computeRMSE(x,batchN, mu0)})


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
    filename = paste("sims/result_",paste(K,balanced ,sigma , initial_truth,sep =  "_"),".Rda",sep="")
    save(result, file=filename)
  }
  # result
}
