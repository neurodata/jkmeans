setwd("~/git/jkmeans/test/")
# build('jkmeans')
# install.packages("jkmeans_1.0.tar.gz", repos = NULL, type = "source")

require("jkmeans")
require(ggplot2)


RMSE <- function(a, b) {
  sqrt(mean((a - b) ^ 2))
}
# MAD <- function(a, b) {
#   median(abs(a - b))
# }


testJK<- function(n, j,   K = 3 , p=1, sigma2 =1) {

  mu0 <- c(1,2,6)
  
  mu <- matrix(rep(mu0, each=n * p), n * K, p)
  
  y <- matrix(rnorm(n * K * p, mu, sigma2), n * K, p)
  
  jk <- jkmeans::jkmeansEM(y, K, j ,  1000, tol = 1E-8)
  
  misclass_error_mat<- matrix(1, n*K, K)
  for(k in 1:K){
    start<- (k-1)*n+1
    end<- k*n
    misclass_error_mat[start:end,k]<- 0
  }

  #deal with label switching
  reOrder<- order(jk$mu[,1])
  mu_est<- jk$mu[ reOrder,]
  
  rmse<- RMSE( mu_est, rep(mu0, p))
  
  mcError<-  sum(misclass_error_mat*jk$zeta[,reOrder])/n/K
  
  c(rmse, mcError)
}

n_seq<- seq(10, 500, by=10)


expiriment <- function(j, m) {
  E3 <- numeric()
  sd_list <- numeric()
  q025<- numeric()
  q975<- numeric()
  
  for (n in n_seq) {
    e <- numeric()
    for (k in 1:m) {
      e <- rbind(e, testJK(n, j, K=3, p = 1 ))
    }
    E3 <- rbind(E3, colMeans(e))
    sd_list<- rbind( sd_list, apply(e, 2,sd))
    q025<- rbind( q025, apply(e, 2, function(x){quantile(x, 0.025)}))
    q975<- rbind( q975, apply(e, 2,function(x){quantile(x, 0.975)}))
    
  }
  list("meanError" = E3, "sdError" = sd_list, "q025Error" = q025, "q0975Error" = q975)
}

E3<- expiriment(j=3,m=100)
E2<- expiriment(j=2,m=100)
E1<- expiriment(j=1,m=100)
# 
# plot(n_seq, E3$meanError[,1],type="l", col="blue")
# lines(n_seq, E3$q025Error[,1], type="l",lty=2, col="blue")
# lines(n_seq, E3$q0975Error[,1], type="l",lty=2, col="blue")
# 
# lines(n_seq, E2$meanError[,1],type="l", col="green")
# lines(n_seq, E2$q025Error[,1], type="l",lty=2, col="green")
# lines(n_seq, E2$q0975Error[,1], type="l",lty=2, col="green")
# 
# lines(n_seq, E1$meanError[,1],type="l", col="red")
# lines(n_seq, E1$q025Error[,1], type="l",lty=2, col="red")
# lines(n_seq, E1$q0975Error[,1], type="l",lty=2, col="red")


plot1<- data.frame( "n" = rep(n_seq, 3),
            "RMSE_avg"=c(E1$meanError[,1], E2$meanError[,1], E3$meanError[,1]),
            "RMSE_2.5q"=c(E1$q025Error[,1], E2$q025Error[,1], E3$q025Error[,1]),
            "RMSE_97.5q"=c(E1$q0975Error[,1], E2$q0975Error[,1], E3$q0975Error[,1]),
            "J"= as.factor(rep(c(1:3),each= length(n_seq)))
            )


p<- ggplot(data=plot1, aes())
p<- p+  geom_line(aes(x= n, y=RMSE_avg, colour=J)) 
p<- p+ geom_line(aes(x= n, y=RMSE_2.5q, colour=J), linetype=2)
p<- p+ geom_line(aes(x= n, y=RMSE_97.5q, colour=J), linetype=2)
p<- p+ theme_bw() + labs(title = "mu = 1,2,6")


pdf("RMSE.pdf",8,6)
p
dev.off()




plot2<- data.frame( "n" = rep(n_seq, 3),
                    "MCerror_avg"=c(E1$meanError[,2], E2$meanError[,2], E3$meanError[,2]),
                    "MCerror_2.5q"=c(E1$q025Error[,2], E2$q025Error[,2], E3$q025Error[,2]),
                    "MCerror_97.5q"=c(E1$q0975Error[,2], E2$q0975Error[,2], E3$q0975Error[,2]),
                    "J"= as.factor(rep(c(1:3),each= length(n_seq)))
)


p<- ggplot(data=plot2, aes())
p<- p+  geom_line(aes(x= n, y=MCerror_avg, colour=J)) 
p<- p+ geom_line(aes(x= n, y=MCerror_2.5q, colour=J), linetype=2)
p<- p+ geom_line(aes(x= n, y=MCerror_97.5q, colour=J), linetype=2)
p<- p+ theme_bw() + labs(title = "mu = 1,2,6")


pdf("ExpectedMisclassification.pdf",8,6)
p
dev.off()



##save a few samples

n<- 100
K<- 3
p<- 1
mu<- c(1,2,6)
sigma2<- 1

y <- matrix(rnorm(n * K * p, rep(mu,each=n*p), sigma2), n * K, p)

#### j=3
jk <- jkmeans::jkmeansEM(y, K, j=3 ,  1000, tol = 1E-8)

reOrder<- order(jk$mu[,1])
mu_est<- jk$mu[ reOrder,]
zeta_est<- jk$zeta[,reOrder]
write.table(y,file="y.csv", row.names = F, col.names  = F, sep = ",")
write.table(mu_est,file="mu_est_j3.csv", row.names = F, col.names  = F, sep = ",")
write.table(zeta_est,file="zeta_est_j3.csv", row.names = F, col.names  = F, sep = ",")




#### j=2
jk <- jkmeans::jkmeansEM(y, K, j=2 ,  1000, tol = 1E-8)

reOrder<- order(jk$mu[,1])
mu_est<- jk$mu[ reOrder,]
zeta_est<- jk$zeta[,reOrder]
write.table(mu_est,file="mu_est_j2.csv", row.names = F, col.names  = F, sep = ",")
write.table(zeta_est,file="zeta_est_j2.csv", row.names = F, col.names  = F, sep = ",")




#### j=1
jk <- jkmeans::jkmeansEM(y, K, j=1 ,  1000, tol = 1E-8)

reOrder<- order(jk$mu[,1])
mu_est<- jk$mu[ reOrder,]
zeta_est<- jk$zeta[,reOrder]
write.table(mu_est,file="mu_est_j1.csv", row.names = F, col.names  = F, sep = ",")
write.table(zeta_est,file="zeta_est_j1.csv", row.names = F, col.names  = F, sep = ",")


