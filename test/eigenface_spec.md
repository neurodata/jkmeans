``` r
require("pixmap")
```

    ## Loading required package: pixmap

``` r
require("jkmeans")
```

    ## Loading required package: jkmeans

``` r
require("mclust")
```

    ## Loading required package: mclust

    ## Package 'mclust' version 5.2

    ## Type 'citation("mclust")' for citing this R package in publications.

``` r
setwd("~/git/jkmeans/lfwcrop_face/lfwcrop_grey/pick/")

faces<- sapply(list.files("."),function(f){
  x=read.pnm(file = f)
  c(x@grey)
})
```

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

    ## Warning in rep(cellres, length = 2): 'x' is NULL so the result will be NULL

``` r
M<- sapply(list.files("."), function(x){substr(x,1,5)})

X<- t(faces)

nf<- 64

K<- 2

plotFace<- function(x){
  image(matrix(x,64,64),col = grey(seq(0, 1, length = 256)))
}

# True Label Count
table(M)
```

    ## M
    ## Alvar Jenni 
    ##    21    11

``` r
gmm<- jkmeansEM(X,k = K,useKmeansIni = T,meansIni = matrix(rep(1,2*10),nrow = 2),sigma2_ini = 0.1,fixW = F,
                tau = 1)

kmeans<- jkmeansEM(X,k = K,useKmeansIni = T,meansIni = matrix(rep(1,2*10),nrow = 2),sigma2_ini = 0.1,fixW = T,
                tau = 1000)

jk<- jkmeansEM(X,k = K,useKmeansIni = T,meansIni = matrix(rep(1,2*10),nrow = 2),sigma2_ini = 0.1,fixW = F,
                tau = 1000)


adjustedRandIndex(gmm$M,M)
```

    ## [1] -0.03254346

``` r
adjustedRandIndex(kmeans$M,M)
```

    ## [1] -0.03254346

``` r
adjustedRandIndex(jk$M,M)
```

    ## [1] -0.03254346
