setwd("~/git/")
# build('jkmeans')
# install.packages("jkmeans_1.0.tar.gz", repos = NULL, type = "source")

require("jkmeans")

n <- 100
p <- 2
K <- 5
sigma2 <- 1

mu <- c(1:K)

mu_post <- c(rnorm(n, mu[1], sqrt(sigma2 / n)),
             rnorm(n, mu[2], sqrt(sigma2 / n)),
             rnorm(n, mu[3], sqrt(sigma2 / n)),
             rnorm(n, mu[4], sqrt(sigma2 / n)))

hist(mu_post, breaks = 100)

y <- matrix(rnorm(n * K, rep(mu, each = n), sigma2))


hist(y, breaks = 100)
# y <- cbind(y,rnorm(n*K),rnorm(n*K))


jk <- jkmeans::jkmeansQNEM(y, K, j = 3,  300)

jk2 <- jkmeans::jkmeansEM(y, K, j = 3,  100)

jk$mu
jk2$mu



# print(jk$mu)
# # print(jk$w)
# print(jk$sigma2)
# hist(jk$zeta)

# image(t(jk$zeta))

