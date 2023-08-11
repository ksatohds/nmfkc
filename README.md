
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
#----------------------------
# obervation and covariates
#----------------------------
library(fda)
data(CanadianWeather)
d <- CanadianWeather$dailyAv[,,1]
Y <- d-min(d)
u0 <- t(CanadianWeather$coordinates[,2:1])
u0[1,] = -u0[1,]
umin <- apply(u0,1,min)
umax <- apply(u0,1,max)
U <- (u0-umin)/(umax-umin) # normalization

#----------------------------
# without covariate
#----------------------------
A <- diag(ncol(Y))
library(nmfreg)
result <- nmfreg(Y,A,Q=2) # Y~XCA=XB

# visualization of some results
plot(result$errs) # check convergence
result$r.squared # coefficient of determination

# check individual fit
n <- 1
plot(Y[,n]) # observation
lines(result$YHAT[,n]) # fitted values

# check basis function of which sum is 1
q <- 1
plot(result$X[,q])

# soft clulustering based on P
Q <- nrow(result$P)
plot(t(U),type="n")
legend("topright",legend=1:Q,fill=1:Q+1)
stars(t(result$P),
      locations=t(U),
      draw.segments = TRUE,labels=NULL,
      col.segments=1:nrow(result$P)+1,
      len=t(U)/10,add=T)

#----------------------------
# with ordinary covariates
#----------------------------
result <- nmfreg(Y,U,Q=2)
result$r.squared # bad fit

#----------------------------
# with covariates using kernel
#----------------------------
# find best seed for cv
best.seed  <- nmfreg.cv.seed(n=ncol(Y),div=10)

# perform cv for some beta
betas <- c(0.1,0.2,0.5,1,2,5,10)
errs <- 0*betas
for(i in 1:length(betas)){
  print(i)
  A <- create.kernel(U,beta=betas[i])
  result <- nmfreg.cv(Y,A,Q=2,div=10,seed=best.seed)
  errs[i] <- result$err
}
# check objective function by beta
plot(betas,errs,type="o",log="x")
(best.beta <- betas[which.min(errs)])

# create kernel with best beta
A <- create.kernel(U,beta=best.beta)
result <- nmfreg(Y,A,Q=2)
result$r.squared # less than nmf without covariates

# soft clulustering based on P
library(akima)
q <- 2
result.interp <- interp(U[1,],U[2,],result$P[q,])
filled.contour(result.interp,xlab=rownames(U)[1],ylab=rownames(U)[2])
```
