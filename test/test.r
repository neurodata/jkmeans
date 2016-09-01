require("jkmeans")

n<- 100

p<- 2
K<- 3
sigma2 <- 1

y<- matrix(0, n*K, p)

mu<- matrix(rep(c(1,2,3),2),K,p)
for (k in 1:K) {
    # print(mu)

    for(i in 1:n){
        idx<- (k-1)*n + i
        y[idx,]<- rnorm(p, mu[k,], sqrt(sigma2))
    }
}


mu_ini<- matrix(rnorm(6),K,p)

jk <- jkmeans::jkmeansEM(y, k=3, j = 3,  1000,tol = 1E-10,useKmeansIni = F, meansIni = mu_ini)


print(jk$mu)
print(jk$M)



# # print(jk$w)
# print(jk$sigma2)
# hist(jk$zeta)

# image(t(jk$zeta))

