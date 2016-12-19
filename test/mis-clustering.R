mean(1/rgamma(10000,2,1))
hist(1/rgamma(10000,2,1))



mrate<-function(y){
test<-sapply(c(1:1000),function(x){
y<- c(y, rnorm(1000,0,0.5))

mu1<- c(1,rnorm(1000,0,0.1))
mu2<- c(2,rnorm(1000,0,0.11))
exp( - (sum((y-mu1)^2) - sum((y-mu2)^2)))
})
sum(test>1)/1000
}

result<-sapply(seq(1.5,5.0,by = 0.1), mrate)

pdf("misrate.pdf",4,3)
plot(seq(1.5,5.0,by = 0.1),result,type = "l")
abline(v=2,lty=2,col="red")
dev.off()
