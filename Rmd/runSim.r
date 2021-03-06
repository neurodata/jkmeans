require('jkmeans')

setwd("~/git/jkmeans/Rmd/")


source("simulation.r")

opt<- rbind(c(2,5), c(T,F), c(0.5,1), c(T,F), c(T,F))



for(i in 1:2){
  for(j in 1:2){
    for(k in 1:2){
      for(l in 1:2){
        simulation(K=opt[1,i],balanced = opt[2,j] ,sigma = opt[3,k], initial_truth= opt[4,l],fix_sigma2 = T)
      }
    }
  }
}

