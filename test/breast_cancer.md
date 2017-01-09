I``` r
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
setwd("~/git/jkmeans/wdbc/")

wdbc<- read.csv("./wdbc.data.txt",header = F,stringsAsFactors = F)

M<- wdbc$V2
X<- wdbc[,3:ncol(wdbc)]
X<- as.matrix(X)

table(M)
```

    ## M
    ##   B   M 
    ## 357 212

``` r
K<- 2

gmm<- jkmeansEM(X,k = K,useKmeansIni = T,meansIni = matrix(rep(1,2*10),nrow = 2),sigma2_ini = 1,fixW = F, tau = 1)

kmeans<- jkmeansEM(X,k = K,useKmeansIni = T,meansIni = matrix(rep(1,2*10),nrow = 2),sigma2_ini = 1,fixW = T, tau = 1000)

jk<- jkmeansEM(X,k = K,useKmeansIni = T,meansIni = matrix(rep(1,2*10),nrow = 2),sigma2_ini = 1,fixW = F,   tau = 1000)

jk_point1<- jkmeansEM(X,k = K,useKmeansIni = T,meansIni = matrix(rep(1,2*10),nrow = 2),sigma2_ini = 1,fixW = F,   tau = 0.1)


gmm$zeta[2,]
```

    ## [1] 2.026621e-16 1.000000e+00

``` r
kmeans$zeta[2,]
```

    ## [1] 0 1

``` r
jk$zeta[2,]
```

    ## [1] 0 1

``` r
jk_point1$zeta[2,]
```

    ## [1] 0.124038 0.875962

``` r
adjustedRandIndex(gmm$M,M)
```

    ## [1] 0.6751701

``` r
adjustedRandIndex(kmeans$M,M)
```

    ## [1] 0.6870801

``` r
adjustedRandIndex(jk$M,M)
```

    ## [1] 0.6751701

``` r
adjustedRandIndex(jk_point1$M,M)
```

    ## [1] 0.7488219
