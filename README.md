# Installation

You can install the development version of nmfkcreg from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ksatohds/nmfkcreg")
```

# Functions

There are three functions in **nmfkcreg** package.

- **nmfkcreg** function optimizes the model
- **nmfkcreg.cv** function is used for k-fold cross-validation
- **create.kernel** function is used for creating kernel matrix from covariates

# Statistical model

1. An ordinary Linear Model (LM) can be written as $y \approx Xb$, where $y$ is the observation, $X=(x_1,...,x_Q)$ is the design matrix, and $b$ is a vector of regression coefficients. 

2. Using the common $X$ among $N$ individual observations, $[y_1,y_2,...,y_N] \approx X[b_1,...,b_N]$，i.e. $Y \approx XB$, this is called Non-Negative Matrix Factorization (NMF). The most significant difference between NMF and LM is that $X$ on NMF is optimized as well as $B$.

3. Since all the matrices are non-negative, the components of the regression coefficients $b_n$ for each $n$ are also non-negative. Therefore we know the contribution of the bases $x_1,.., x_Q$ and those probabilities $p_{1,n},...,p_{Q,n}$ can be used for clustering.

4. Furthermore, the regression coefficient matrix can be explained by covariates and its Gaussian kernel can be applied in Satoh (2023).
Formally the model is contained in tri-NMF by Ding et al. (1964) and the update formula for optimization has been derived. The model can be described as a product of three matrices and is related to the growth curve model by Potthoff and Roy (1964). 

# Matrices

The goal of **nmfkcreg** is to optimize $X(P,Q)$ and $C(Q,N)$ on the NMF (Non-negative Matrix Factorization) with kernel covariates regression, $Y(P,N) \approx X(P,Q)B(Q,N)=X(P,Q)C(Q,N)A(N,N)$ where $Y(P,N)$ and $A(N,N)$ are given.

- $Y(P,N)=(y_1,...y_N)$: **given** observation matrix
- $A(N,N)$: **given** kernel matrix of which $(i,j)$ element can be 
written as $K(u_i,u_j)=exp(−\beta|u_i-u_j|^2)$ here $U(R,N)=(u_1,...u_N)$ 
is covariate matrix. Note that identity matrix is used when there are no covariates. Or matrix $A(R,N)$ having N columns can be accepted.
- $X(P,Q)$: **unknown** basis matrix whose column sum is 1 and Q<=min{P,N}.
Q is the number of basis (rank).
-  $C(Q,N)$: **unknown** parameter matrix which is described by $\Theta$ in the paper Satoh (2023). 
- $B(Q,N)=C(Q,N)A(N,N)$ is regression coefficient matrix.
- $P(Q,N)$ is probability matrix whose column sum is 1 for soft clustering based on regression coefficient matrix $B(Q,N)$.

# References

- Kenichi Satoh (2023)
  On Non-negative Matrix Factorization Using Gaussian Kernels as Covariates,
  Japanese Journal of Applied Statistics, 52 (2), in press.
[(preprint)](https://drive.google.com/file/d/1MnbJOPlcm0hn27WpP8rvcAgzy5X2E53B/view?usp=sharing)
- Ding, Chris and Li, Tao and Peng, Wei and Park, Haesun (2006)
  Orthogonal Nonnegative Matrix Tri-Factorizations for Clustering,
  Proceedings of the 12th ACM SIGKDD international conference on Knowledge discovery and data mining, 126-135.
  https://doi.org/10.1145/1150402.1150420
- Potthoff, Richard F., and Roy, S. N. (1964) 
  A generalized multivariate analysis of variance model useful especially for growth curve problems,
  Biometrika, 51 (3/4), 313–326.
  https://doi.org/10.2307/2334137

# Examples

1. Longitudinal data: COVID-19 in Japan
2. Topic model: data_corpus_inaugural
3. Spatiotemporal Analysis: CanadianWeather
4. Kernel ridge regression: mcycle

## 1. Longitudinal data
- COVID-19 in Japan
- https://www3.nhk.or.jp/news/special/coronavirus/data/
``` r
library(nmfkcreg)
d <- read.csv(
  "https://www3.nhk.or.jp/n-data/opendata/coronavirus/nhk_news_covid19_prefectures_daily_data.csv")
colnames(d) <- c(
  "Date","Prefecture_code","Prefecture_name",
  "Number_of_infected","Cumulative_Number_of_infected",
  "Number_of_deaths","Cumulative_Number_of_deaths",
  "Number_of_infected_100000_population_in_the_last_week")
n <- length(unique(d$Prefecture_code)) # 47
Y <- matrix(d$Number_of_infected,nrow=nrow(d)/n,ncol=n)
colnames(Y) <- unique(d$Prefecture_code)
rownames(Y) <- unique(d$Date)
Y <- Y[rowSums(Y)!=0,]

# nmf
Q <- 7
result <- nmfkcreg(Y,Q=Q)
result$r.squared # goodness of fit

# individual fit
n <- 1
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1,cex=1)
plot(Y[,n],type="l",main=n)
lines(result$XB[,n],col=2)

# basis function of which sum is 1
par(mfrow=c(Q,1),mar=c(0,0,0,0),cex=1)
for(q in 1:Q){
  barplot(result$X[,q],col=q+1,border=q+1,las=3,
    ylim=range(result$X),ylab=paste0("topic",q)) 
  legend("left",fill=q+1,legend=q)
}

# individual coefficients
n <- 1
result$B[,n]

# individual probabilites on basis function
n <- 1
result$B[,n]/sum(result$B[,n])
result$P[,n]

# hard clustering based on P
mycol <- result$cluster
library(NipponMap)
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
JapanPrefMap(col=mycol+1)
legend("left",fill=1:Q+1,legend=1:Q)
``` 

## 2. Topic model
- data_corpus_inaugural
- US presidential inaugural address texts
``` r
library(nmfkcreg)
#------------------
# text analysis
#------------------
library(quanteda)
corp <- corpus(data_corpus_inaugural)
tok <- tokens(corp)
tok <- tokens_remove(tok,pattern=stopwords("en",source="snowball"))
df <- dfm(tok)
df <- dfm_select(df,min_nchar=3)
df <- dfm_trim(df,min_termfreq=100)
d <- as.matrix(df)
index <- order(colSums(d),decreasing=T) 
d <- d[,index] # document-word matrix
colSums(d)[1:30] # Top 30 most frequent words

#------------------
# without covariates
#------------------
Y <- t(d)
Q <- 3
result <- nmfkcreg(Y,Q=Q)
result$r.squared # coefficient of determination
colnames(result$P) <- corp$Year
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1,cex=1)
barplot(result$P,col=1:Q+1,legend=T,las=3,
  ylab="Probabilities of topics")

# basis function of which sum is 1
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=0.6)
barplot(t(result$X),col=1:Q+1,las=3,
  ylab="Probabilities of words on each topic") 
legend("topright",fill=1:Q+1,legend=paste0("topic",1:Q))

# contribution of words to each topics
Xp <- prop.table(result$X,1)
head(rowSums(Xp))
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=0.6)
barplot(t(Xp),las=3,col=1:Q+1,
  ylab="Proportion of words on each topic")
legend("topright",fill=1:Q+1,legend=paste0("topic",1:Q))
Xp[Xp[,1]>0.7,] # featured words on topic1
Xp[Xp[,2]>0.7,] # featured words on topic2
Xp[Xp[,3]>0.7,] # featured words on topic3

#------------------
# with covariates using covariate matrix U
#------------------
U <- t(as.matrix(corp$Year))
# k-fold cross validation for beta
betas <- c(0.2,0.5,1,2,5)/10000
objfuncs <- 0*betas
for(i in 1:length(betas)){
  print(i)
  A <- create.kernel(U,beta=betas[i])
  result <- nmfkcreg.cv(Y,A,Q,div=5)
  objfuncs[i] <- result$objfunc
}
table(result$block) # partition block of cv
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(betas,objfuncs,type="o",log="x")
(best.beta <- betas[which.min(objfuncs)])

# create kernel with best beta
A <- create.kernel(U,beta=best.beta)
result <- nmfkcreg(Y,A,Q)
result$r.squared # less than nmf without covariates

# Topic probability changing over time
colnames(result$P) <- corp$Year
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
barplot(result$P,col=1:Q+1,legend=T,las=3,ylab="Probability of topic")
``` 

## 3. Spatiotemporal Analysis
-  CanadianWeather
``` r
library(nmfkcreg)
library(fda)
data(CanadianWeather)
d <- CanadianWeather$dailyAv[,,1]
Y <- d-min(d)
u0 <- CanadianWeather$coordinates[,2:1]
u0[,1] <- -u0[,1]

#------------------
# without covariate
#------------------
result <- nmfkcreg(Y,Q=2)
result$r.squared # coefficient of determination

# visualization of some results
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(result$objfunc.iter) # convergence

# individual fit
par(mfrow=c(6,6),mar=c(0,0,0,0)+0.1,cex=1)
for(n in 1:ncol(Y)){
  plot(Y[,n],ylim=range(Y),axes=F,type="l") # observation
  lines(result$XB[,n],col=2) # fitted values
  abline(h=-min(d),col="gray",lty=3,lwd=3)
  text(median(1:365),0,colnames(Y)[n],pos=3)
  box()
}

# basis function of which sum is 1
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(result$X[,1],type="n",ylim=range(result$X[,1]),
  ylab="basis function")
Q <- ncol(result$X)  
for(q in 1:Q) lines(result$X[,q],col=q+1)
legend("topright",legend=1:Q,fill=1:Q+1)

# soft clustering based on P
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(u0,type="n")
legend("topright",legend=1:Q,fill=1:Q+1)
stars(t(result$P),
      locations=u0,scale=F,
      draw.segments=TRUE,labels=colnames(Y),
      col.segments=1:Q+1,
      len=max(u0)/30,add=T)

#------------------
# with covariates using location information
#------------------
u <- t(u0)
umin <- apply(u,1,min)
umax <- apply(u,1,max)
U <- (u-umin)/(umax-umin) # normalization
result <- nmfkcreg(Y,A=U,Q=2)
result$r.squared

#------------------
# with covariates using kernel matrix A
#------------------
# U=[u1,...,uN]
# A= |K(u1,u1),...,K(u1,uN)|
#    |                     |
#    |K(uN,u1),...,K(uN,uN)|
# K(u,v)=exp{-beta*|u-v|^2}

# k-fold cross validation for beta
betas <- c(0.5,1,2,5,10)
objfuncs <- 0*betas
for(i in 1:length(betas)){
  print(i)
  A <- create.kernel(U,beta=betas[i])
  result <- nmfkcreg.cv(Y,A,Q=2,div=10)
  objfuncs[i] <- result$objfunc
}
table(result$block) # partition block of cv

# objective function by beta
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(betas,objfuncs,type="o",log="x")
(best.beta <- betas[which.min(objfuncs)])

# create kernel with best beta
A <- create.kernel(U,beta=best.beta)
result <- nmfkcreg(Y,A,Q=2)
result$r.squared # less than nmf without covariates

# prediction of coefficients(b) on mesh point V
# U=[u1,...,uN]
# V=[v1,...,vM]
# A= |K(u1,v1),...,K(u1,vM)|
#    |                     |
#    |K(uN,v1),...,K(uN,vM)|
v <- seq(from=0,to=1,length=20)
V <- t(cbind(expand.grid(v,v)))
dim(V)
A <- create.kernel(U,V,beta=best.beta)
B <- result$C %*% A
P <- prop.table(B,2)
P[,1:6]

# soft clustering based on P (base func 2) by using covariates
z <- matrix(P[2,],nrow=length(v)) 
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
filled.contour(v,v,z,
               plot.axes={
                 points(t(U),col=3,pch=19)
                 text(t(U),colnames(U),pos=4)
               }
)
```

## 4. Kernel ridge regression
- mcycle
``` r
library(nmfkcreg)
library(MASS)
d <- mcycle
x <- d$times
y <- d$accel
Y <- t(as.matrix(y-min(y)))
dim(Y)
U <- t(as.matrix(x))
dim(U)
# scatter plot
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(U,Y)

# linear curve
A <- rbind(1,U)
dim(A)
result <- nmfkcreg(Y,A,Q=1)
result$r.squared
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(U,Y)
lines(U,result$XB,col=2)

# cv for optimization of beta and gamma
betas <- c(1,2,5,10,20)/100
objfuncs <- 0*(1:length(betas))
for(i in 1:length(betas)){
  print(i)
  A <- create.kernel(U,beta=betas[i])
  result <- nmfkcreg.cv(Y,A,Q=1,div=10)
  objfuncs[i] <- result$objfunc
}
# objective function by beta
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(betas,objfuncs,type="o",log="x")
table(result$block) # partition block of cv

(beta.best <- betas[which.min(objfuncs)])  
A <- create.kernel(U,beta=beta.best)
result <- nmfkcreg(Y,A,Q=1)
result$r.squared

# fitted curve
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(x,as.vector(Y))
lines(x,as.vector(result$XB),col=2)
```

# Author
-  Kenichi Satoh, [homepage](https://sites.google.com/view/ksatoh/english)
