
# nmfreg

<!-- badges: start -->
<!-- badges: end -->

The goal of nmfreg is to perform NMF regression model described by Y~XCA where 
  observation matrix Y:PxN,
  kernel matrix A:NxN with parameter beta,
  rank of basis matrix X:PxQ whose column sum is 1
  and coefficient matrix C:QxN,
  Y and A are known, and X and C are unknown.

## Installation

You can install the development version of nmfreg from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ksatohds/nmfreg")
```

## Example

``` r
#----------------------------
# example 1: iris
## without covariate
#----------------------------
Y <- t(iris[,-5])
A <- diag(ncol(Y))
library(nmfreg)
result <- nmfreg(Y,A,Q=2) # Y~XCA=XB

# visualization of some results
plot(result$objfunc.iter) # check convergence

# check fit
plot(as.vector(result$YHAT),as.vector(result$Y),
main=paste0("r.squared=",round(result$r.squared,3)))
abline(a=0,b=1,col=2)

# soft clulustering based on P
labels <- as.numeric(iris[,5])
plot(result$P[1,],col=labels)
legend("topright",
  legend=unique(iris[,5]),fill=unique(labels))

#----------------------------
# example 2: basketball players and statistics
## without covariate
#----------------------------
# https://rpubs.com/sirinya/847402
d <- read.csv("http://datasets.flowingdata.com/ppg2008.csv")
y <- d[,-1]
rownames(y) <- d[,1]

Y <- t(y)
A <- diag(ncol(Y))
library(nmfreg)
result <- nmfreg(Y,A,Q=2) # Y~XCA=XB

# visualization of some results
plot(result$objfunc.iter) # check convergence

# check fit
plot(as.vector(result$YHAT),as.vector(result$Y),
     main=paste0("r.squared=",round(result$r.squared,3)))
abline(a=0,b=1,col=2)

# check individual fit
n <- 1
f <- rbind(Y[,n],result$YHAT[,n])
rownames(f) <- c("obs","fitted")
barplot(f,beside=T,las=3,legend=T) # observation

# check basis function of which sum is 1
Q <- ncol(result$X)
barplot(t(result$X),beside=T,col=1:Q+1,legend=T,las=3)

# soft clulustering based on P
stars(t(result$P),scale=F,
      draw.segments=TRUE,labels=colnames(Y),
      col.segments=1:Q+1,
      len=1)

#----------------------------
# example 3: CanadianWeather
## 3.1 without covariate
## 3.2 with covariates
## 3.3 with covariates using kernel
#----------------------------
library(fda)
data(CanadianWeather)
d <- CanadianWeather$dailyAv[,,1]
Y <- d-min(d)
u0 <- CanadianWeather$coordinates[,2:1]
u0[,1] = -u0[,1]
u = t(u0)
umin <- apply(u,1,min)
umax <- apply(u,1,max)
U <- (u-umin)/(umax-umin) # normalization
library(nmfreg)

#------------------
## 3.1 without covariate
#------------------
A <- diag(ncol(Y))
result <- nmfreg(Y,A,Q=2) # Y~XCA=XB

# visualization of some results
plot(result$objfunc.iter) # check convergence
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
plot(u0,type="n")
legend("topright",legend=1:Q,fill=1:Q+1)
stars(t(result$P),
      locations=u0,scale=F,
      draw.segments=TRUE,labels=colnames(Y),
      col.segments=1:Q+1,
      len=max(u0)/30,add=T)

#------------------
## 3.2 with covariates
#------------------
result <- nmfreg(Y,U,Q=2) # Y~XCA=XB
result$r.squared # bad fit

#------------------
## 3.3 with covariates using kernel
#------------------
# perform cv for some beta
betas <- c(0.5,1,2,5,10)
objfuncs <- 0*betas
for(i in 1:length(betas)){
  print(i)
  A <- create.kernel(U,beta=betas[i])
  result <- nmfreg.cv(Y,A,Q=2,div=10)
  objfuncs[i] <- result$objfunc
}
table(result$block) # partition block of cv

# check objective function by beta
plot(betas,objfuncs,type="o",log="x")
(best.beta <- betas[which.min(objfuncs)])

# create kernel with best beta
A <- create.kernel(U,beta=best.beta)
result <- nmfreg(Y,A,Q=2) # Y~XCA=XB
result$r.squared # less than nmf without covariates

# soft clulustering based on P
library(akima)
q <- 2
result.interp <- interp(U[1,],U[2,],result$P[q,])
filled.contour(result.interp,
  xlab=rownames(U)[1],ylab=rownames(U)[2])
```
