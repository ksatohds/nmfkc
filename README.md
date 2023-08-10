
# nmfreg

<!-- badges: start -->
<!-- badges: end -->

The goal of nmfreg is to ...

## Installation

You can install the development version of nmfreg from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ksatohds/nmfreg")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(fda)
data(CanadianWeather)
d <- CanadianWeather$dailyAv[,,1]
Y <- d-min(d)

# without covariate
A <- diag(ncol(Y))
library(nmfreg)
result <- nmfreg(Y,A,Q=2)
result$r.squared
n <- 1
plot(Y[,n])
lines(result$YHAT[,n])

# with covariates
u0 <- t(CanadianWeather$coordinates[,2:1])
u0[1,] = -u0[1,]
umin <- apply(u0,1,min)
umax <- apply(u0,1,max)
U <- (u0-umin)/(umax-umin) # normalization

myseed  <- nmfreg.cv.seed(n=ncol(Y),div=10)
betas <- c(0.1,0.2,0.5,1,2,5,10)
errs <- 0*betas
for(i in 1:length(betas)){
  print(i)
  A <- create.kernel(U,beta=betas[i])
  result <- nmfreg.cv(Y,A,Q=2,div=10,seed=myseed)
  errs[i] <- result$err
}
plot(betas,errs,type="o",log="x")
```

