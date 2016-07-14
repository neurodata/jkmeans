setwd("~/git/")
# build('jkmeans')
# install.packages("jkmeans_1.0.tar.gz", repos = NULL, type = "source")

require("jkmeans")
require(ggplot2)


RMSE <- function(a, b) {
  sqrt(mean((a - b) ^ 2))
}
MAD <- function(a, b) {
  median(abs(a - b))
}


testJK<- function(n, j,   K = 3 , p=3, sigma2 =1) {

  mu0 <- c(1,2,3)
  
  mu <- matrix(rep(c(1:K), n * p), n * K, p)
  
  y <- matrix(rnorm(n * K * p, mu, sigma2), n * K, p)
  
  jk <- jkmeans::jkmeansEM(y, K, j ,  1000, tol = 1E-8)
  
  mu_est<- jk$mu[ order(jk$mu[,1]),]
  
  c(RMSE( mu_est, rep(mu0, p)), MAD(mu_est , rep(mu0, p)))
  
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
      e <- rbind(e, testJK(n, j ))
    }
    E3 <- rbind(E3, colMeans(e))
    sd_list<- rbind( sd_list, apply(e, 2,sd))
    q025<- rbind( q025, apply(e, 2, function(x){quantile(x, 0.025)}))
    q975<- rbind( q975, apply(e, 2,function(x){quantile(x, 0.975)}))
    
  }
  list("meanError" = E3, "sdError" = sd_list, "q025Error" = q025, "q0975Error" = q975)
}

E3<- expiriment(3,100)
E2<- expiriment(2,100)
E1<- expiriment(1,100)
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
p<- p+ theme_bw() + labs(title = "mu = (1,1,1), (2,2,2), (3,3,3)")


pdf("test2.pdf",8,6)
p
dev.off()
