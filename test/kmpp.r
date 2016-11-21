kmpp <- function(X, k) {
  n <- nrow(X)
  C <- numeric(k)
  C[1] <- sample(1:n, 1)
  
  for (i in 2:k) {
    dm <- abs(outer(X, X[C, ], "-"))
    pr <- apply(dm, 1, min)
    pr[C] <- 0
    C[i] <- sample(1:n, 1, prob = pr)
  }
  
  X[C, ]
}