kmpp <- function(X, k) {
  n <- nrow(X)
  C <- numeric(k)
  C[1] <- sample(1:n, 1)
  
  D2<-  apply(X,1,function(x){       sum((x-X[C[1],])^2)})
  
  for (i in 2:k) {
    
    pr<- D2
    pr[C] <- 0
    C[i] <- sample(1:n, 1, prob = pr)
    
    D2new <-  apply(X,1,function(x){       sum((x-X[C[i],])^2)})
    D2<- pmin(D2new,D2)
  }
  
  X[C, ]
}