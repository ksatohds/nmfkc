
# nmfreg

<!-- badges: start -->
<!-- badges: end -->

The goal of **nmfreg** is to optimize X and C on the NMF regression model $Y \approx X C A$.

- given Y[P,N] observation matrix
- given A[N,N] kernel matrix of which $(i,j)$ element can be 
written as $K(u_i,u_j)=exp(−\beta|u_i-u_j|^2)$ here $U=[u_1,...u_N]$ 
is covariate matrix. Note that identity matrix is used as A=diag(ncol(Y)) in the case where there are no covariate.
- unknown X[P,Q] basis matrix whose column sum is 1.
Q shows the number of basis (rank).
- unkown C[Q,N] regression coefficient matrix which is described by $\Theta$ in the paper Satoh (2023).
- **nmfreg** function is used for optimization of X and C
- **nmfreg.cv** function is used for k-fold cross-validation optimizing $\beta$ and $Q$
- **create.kernel** function is used fot creating kernel matrix A from covariate matrix U
## Reference

- Satoh, K. (2023) On Non-negative Matrix Factorization Using Gaussian Kernels as Covariates, Japanese Journal of Applied Statistics 52 (2), in press. (in Japanese)
- Satoh, K. (2022) Soft Clustering Based on Non-negative Matrix
Factorization for Longitudinal Data, {\it Japanese Journal of Applied Statistics}, {\bf 51} (1&2), 1-18. https://doi.org/10.5023/jappstat.51.1 (in Japanese)
- Ding, C., Tao, L., Wei, P. and Haesun, P. (2006)
Orthogonal Nonnegative Matrix Tri-Factorizations for Clustering,
 {\it Proceedings of the 12th ACM SIGKDD international conference on Knowledge discovery and data mining}, 126-135.
- Potthoff, R. F. and Roy, S. N. (1964) A generalized multivariate analysis of variance model useful especially for growth curve problems,
{\it Biometrika}, {\bf 51}, 313–326.


## Reference (in Japanese)
- 佐藤健一 (2023) ガウスカーネルを共変量に用いた非負値行列因子分解について, 応用統計学 52 (2), 印刷中.
- 佐藤健一 (2022) 経時測定データに対する非負値行列因子分解によるソフトクラスタリングについて, 応用統計学, 51(1-2), 1-18. https://doi.org/10.5023/jappstat.51.1

## Installation

You can install the development version of nmfreg from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ksatohds/nmfreg")
```

## Example

There are three examples.

### iris
``` r
library(nmfreg)
Y <- t(iris[,-5])
A <- diag(ncol(Y))
result <- nmfreg(Y,A,Q=2) # Y~XCA=XB

# visualization of some results
plot(result$objfunc.iter) # check convergence

# check fit
plot(as.vector(result$YHAT),as.vector(result$Y),
main=paste0("r.squared=",round(result$r.squared,3)))
abline(a=0,b=1,col=2)

# dimension reduction based on regression coefficient B
labels <- as.numeric(iris[,5])
plot(t(result$B),col=labels)
legend("topright",
  legend=unique(iris[,5]),fill=unique(labels))
``` 

### basketball players and statistics
- https://rpubs.com/sirinya/847402
``` r
library(nmfreg)
d <- read.csv("http://datasets.flowingdata.com/ppg2008.csv")
y <- d[,-1]
rownames(y) <- d[,1]

Y <- t(y)
A <- diag(ncol(Y))
result <- nmfreg(Y,A,Q=2) # Y~XCA=XB

# visualization of some results
plot(result$objfunc.iter) # check convergence

# good of fit
plot(as.vector(result$YHAT),as.vector(result$Y),
     main=paste0("r.squared=",round(result$r.squared,3)))
abline(a=0,b=1,col=2)

# individual fit
n <- 1
f <- rbind(Y[,n],result$YHAT[,n])
rownames(f) <- c("obs","fitted")
barplot(f,beside=T,las=3,legend=T,main=colnames(Y)[n])

# basis function of which sum is 1
Q <- ncol(result$X)
barplot(t(result$X),beside=T,col=1:Q+1,legend=T,las=3)

# soft clulustering based on P
stars(t(result$P),scale=F,
      draw.segments=TRUE,labels=colnames(Y),
      col.segments=1:Q+1,
      len=1)
``` 

### CanadianWeather
- common preparation
- without covariate
- with covariates using covariate matrix U
- with covariates using kernel matrix A

#### common preparation
``` r
library(nmfreg)
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
``` 

#### without covariate
``` r
A <- diag(ncol(Y))
result <- nmfreg(Y,A,Q=2) # Y~XCA=XB

# visualization of some results
plot(result$objfunc.iter) # check convergence
result$r.squared # coefficient of determination

# individual fit
n <- 1
plot(Y[,n],main=colnames(Y)[n]) # observation
lines(result$YHAT[,n],col=2) # fitted values
legend("topright",
  legend=c("obs","fitted"),fill=c(1,2))

# basis function of which sum is 1
plot(result$X[,1],type="n",ylim=range(result$X[,1]),
  ylab="basis function")
Q <- ncol(result$X)  
for(q in 1:Q) lines(result$X[,q],col=q+1)
legend("topright",legend=1:Q,fill=1:Q+1)

# soft clulustering based on P
plot(u0,type="n")
legend("topright",legend=1:Q,fill=1:Q+1)
stars(t(result$P),
      locations=u0,scale=F,
      draw.segments=TRUE,labels=colnames(Y),
      col.segments=1:Q+1,
      len=max(u0)/30,add=T)
```

#### with covariates using covariate matrix U
``` r
result <- nmfreg(Y,U,Q=2) # Y~XCA=XB
result$r.squared # bad fit
```

#### with covariates using kernel matrix A
``` r
# k-fold cross validation for beta
betas <- c(0.5,1,2,5,10)
objfuncs <- 0*betas
for(i in 1:length(betas)){
  print(i)
  A <- create.kernel(U,beta=betas[i])
  result <- nmfreg.cv(Y,A,Q=2,div=10)
  objfuncs[i] <- result$objfunc
}
table(result$block) # partition block of cv

# objective function by beta
plot(betas,objfuncs,type="o",log="x")
(best.beta <- betas[which.min(objfuncs)])

# create kernel with best beta
A <- create.kernel(U,beta=best.beta)
result <- nmfreg(Y,A,Q=2) # Y~XCA=XB
result$r.squared # less than nmf without covariates

# soft clulustering based on P by using covariates
library(akima)
q <- 2
result.interp <- interp(U[1,],U[2,],result$P[q,])
filled.contour(result.interp,
  xlab=rownames(U)[1],ylab=rownames(U)[2])
```
