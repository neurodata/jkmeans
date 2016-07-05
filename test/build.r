library(devtools)

setwd("~/git/")
build('jkmeans')
install.packages("jkmeans_1.0.tar.gz", repos = NULL, type = "source")

require("jkmeans")

n<- 5
p<- 2
K<- 4
sigma2 <- 1

y<- matrix(0, n*K, p)


for (k in 1:K) {
    mu<- rnorm(p)
    print(mu)

    for(i in 1:n){
        idx<- (k-1)*n + i
        y[idx,]<- rnorm(p, mu, sqrt(sigma2))
    }
}

hist(y, breaks = 30)

# y[,2]<-0




# y

jk<- jkmeans::jkmeans(y, 4, j = 4,  10000)
print(jk)



