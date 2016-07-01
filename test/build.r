library(devtools)

setwd("~/git/")
# build('jkmeans')
# install.packages("jkmeans_1.0.tar.gz", repos = NULL, type = "source")

require("jkmeans")

n<- 10
p<- 10
K<- 3
sigma2 <- 1

y<- matrix(0, n*K, p)

for (k in 1:K) {
    mu<- rnorm(p)*5
    # print(mu)

    for(i in 1:n){
        idx<- (k-1)*n + i
        y[idx,]<- rnorm(p, mu, sqrt(sigma2))
    }
}


# y


jk<- jkmeans::jkmeans(y, K, j = 2,  10000)
jk
