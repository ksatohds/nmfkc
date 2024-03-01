# Installation

You can install the development version of nmfkc from [GitHub](https://github.com/) with:

``` r
#install.packages("remotes")
remotes::install_github("ksatohds/nmfkc")
```

# Functions

There are three functions in **nmfkc** package.

- **nmfkc** function optimizes the model
- **nmfkc.cv** function is used for k-fold cross-validation
- **nmfkc.kernel** function is used for creating kernel matrix from covariates

# Statistical model

1. An ordinary Linear Model (LM) can be written as $y \approx Xb$, where $y$ is the observation, $X=(x_1,...,x_Q)$ is the design matrix, and $b$ is a vector of regression coefficients. 

2. Using the common $X$ among $N$ individual observations, $[y_1,y_2,...,y_N] \approx X[b_1,...,b_N]$，i.e. $Y \approx XB$, this is called Non-Negative Matrix Factorization (NMF). The most significant difference between NMF and LM is that $X$ on NMF is optimized as well as $B$.

3. Since all the matrices are non-negative, the components of the regression coefficient $b=(b_1,...,b_Q)$ are also non-negative. Therefore the proportion $p_q=b_q/(b_1+...+b_Q)$ can be used for clustering.

4. Furthermore, the coefficient matrix can be explained by covariates.
Satoh (2023) proposed Gaussian kernel function as covariates.
Formally the model is contained in tri-NMF by Ding et al. (1964) and the update formula for optimization has been derived. The model can be described as a product of three matrices and is related to the growth curve model by Potthoff and Roy (1964). 

# Matrices

The goal of **nmfkc** is to optimize $X(P,Q)$ and $C(Q,R)$ on the Non-negative Matrix Factorization, $Y(P,N) \approx X(P,Q)C(Q,R)A(R,N)=XB(Q,N)$ where $Y(P,N)$ and $A(R,N)$ are given.

- $Y(P,N)=(y_1,...y_N)$: **given** observation matrix
- $A(R,N)$: **given** covariate matrix. 
  The kernel matrix A(N,N) can be created by nmfkc.kernel function and its $(i,j)$ element
  can be written as Gauss kernel, $K(u_i,u_j)=exp(−\beta|u_i-u_j|^2)$ here $U(R,N)=(u_1,...u_N)$ 
is covariate matrix. Note that identity matrix is used when there are no covariates.
- $X(P,Q)$: **unknown** basis matrix. Q is the number of basis (rank) and Q<=min(P,N).
-  $C(Q,R)$: **unknown** parameter matrix which is described by $\Theta$ in the paper Satoh (2023). 
- $B(Q,N)=C(Q,R)A(R,N)$ is coefficient matrix.

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
2. Spatiotemporal Analysis: CanadianWeather
3. Topic model: data_corpus_inaugural
4. Kernel ridge regression: mcycle
5. Growth curve model: Orthodont

## 1. Longitudinal data
- COVID-19 in Japan
- https://www3.nhk.or.jp/news/special/coronavirus/data/
``` r
library(nmfkc)
d <- read.csv(
  "https://www3.nhk.or.jp/n-data/opendata/coronavirus/nhk_news_covid19_prefectures_daily_data.csv")
colnames(d) <- c(
  "Date","Prefecture_code","Prefecture_name",
  "Number_of_infected","Cumulative_Number_of_infected",
  "Number_of_deaths","Cumulative_Number_of_deaths",
  "Number_of_infected_100000_population_in_the_last_week")
n <- length(unique(d$Prefecture_code)) # 47
Y <- matrix(d$Number_of_infected,nrow=nrow(d)/n,ncol=n)
colnames(Y) <- unique(d$Prefecture_name)
rownames(Y) <- unique(d$Date)
Y <- Y[rowSums(Y)!=0,]

# nmf
Q <- 7
result <- nmfkc(Y,Q=Q)
result$r.squared # goodness of fit

# individual fit
par(mfrow=c(7,7),mar=c(0,0,0,0)+0.1,cex=1)
for(n in 1:ncol(Y)){
  plot(Y[,n],axes=F,type="l") # observation
  lines(result$XB[,n],col=2) # fitted values
  text(1,max(Y[,n])/2,colnames(Y)[n],pos=4)
  box()
}

# basis function of which sum is 1
par(mfrow=c(Q,1),mar=c(0,0,0,0),cex=1)
for(q in 1:Q){
  barplot(result$X[,q],col=q+1,border=q+1,las=3,
    ylim=range(result$X),ylab=paste0("topic",q)) 
  legend("left",fill=q+1,legend=q)
}

# cluster membership probability based on coefficients
n <- 1
result$B[,n]
result$B[,n]/sum(result$B[,n])
result$B.prob[,n]

# hard clustering based on B.prob
mycol <- result$B.cluster
library(NipponMap)
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
JapanPrefMap(col=mycol+1)
legend("left",fill=1:Q+1,legend=1:Q)
``` 


## 2. Spatiotemporal Analysis
-  CanadianWeather
``` r
library(nmfkc)
library(fda)
data(CanadianWeather)
d <- CanadianWeather$dailyAv[,,1]
Y <- d-min(d)
u0 <- CanadianWeather$coordinates[,2:1]
u0[,1] <- -u0[,1]

#------------------
# without covariate
#------------------
result <- nmfkc(Y,Q=2)
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

# cluster membership probability based on coefficients
n <- 1
result$B[,n]
result$B[,n]/sum(result$B[,n])
result$B.prob[,n]

# soft clustering based on B.prob
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(u0,type="n")
legend("topright",legend=1:Q,fill=1:Q+1)
stars(t(result$B.prob),
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
A <- rbind(rep(1,ncol(Y)),U)
result <- nmfkc(Y,A,Q=2)
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
  A <- nmfkc.kernel(U,beta=betas[i])
  result <- nmfkc.cv(Y,A,Q=2,div=10)
  objfuncs[i] <- result$objfunc
}
table(result$block) # partition block of cv

# objective function by beta
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(betas,objfuncs,type="o",log="x")
(best.beta <- betas[which.min(objfuncs)])

# create kernel with best beta
A <- nmfkc.kernel(U,beta=best.beta)
result <- nmfkc(Y,A,Q=2)
result$r.squared # less than nmf without covariates

# prediction of coefficients(b) on mesh point V
# U=[u1,...,uN]
# V=[v1,...,vM]
# A= |K(u1,v1),...,K(u1,vM)|
#    |                     |
#    |K(uN,v1),...,K(uN,vM)|
v <- seq(from=0,to=1,length=20)
V <- t(cbind(expand.grid(v,v)))
plot(t(V))
A <- nmfkc.kernel(U,V,beta=best.beta)
B <- result$C %*% A
P <- prop.table(B,2)
P[,1:6]

# soft clustering based on B.prob (basis function 2) by using covariates
z <- matrix(P[2,],nrow=length(v)) 
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
filled.contour(v,v,z,main="probability of basis function 2",
  color.palette = function(n) hcl.colors(n,"Greens3",rev=TRUE),
               plot.axes={
                 points(t(U),col=7,pch=19)
                 text(t(U),colnames(U),pos=3)
               }
)
```

## 3. Topic model
- data_corpus_inaugural
- US presidential inaugural address texts
``` r
library(nmfkc)
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
result <- nmfkc(Y,Q=Q)
result$r.squared # coefficient of determination
colnames(result$B.prob) <- corp$Year
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1,cex=1)
barplot(result$B.prob,col=1:Q+1,legend=T,las=3,
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
  A <- nmfkc.kernel(U,beta=betas[i])
  result <- nmfkc.cv(Y,A,Q,div=5)
  objfuncs[i] <- result$objfunc
}
table(result$block) # partition block of cv
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(betas,objfuncs,type="o",log="x")
(best.beta <- betas[which.min(objfuncs)])

# create kernel with best beta
A <- nmfkc.kernel(U,beta=best.beta)
result <- nmfkc(Y,A,Q)
result$r.squared # less than nmf without covariates

# Topic probability changing over time
colnames(result$B.prob) <- corp$Year
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
barplot(result$B.prob,col=1:Q+1,legend=T,las=3,ylab="Probability of topic")
``` 


## 4. Kernel ridge regression
- mcycle
``` r
library(nmfkc)
library(MASS)
d <- mcycle
x <- d$times
y <- d$accel
Y <- t(as.matrix(y-min(y)))
U <- t(as.matrix(x))
# scatter plot
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(U,Y)

# linear curve
A <- rbind(1,U)
result <- nmfkc(Y,A,Q=1)
result$r.squared
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(U,Y)
lines(U,result$XB,col=2)

# cv for optimization of beta and gamma
betas <- c(1,2,5,10,20)/100
objfuncs <- 0*(1:length(betas))
for(i in 1:length(betas)){
  print(i)
  A <- nmfkc.kernel(U,beta=betas[i])
  result <- nmfkc.cv(Y,A,Q=1,div=10)
  objfuncs[i] <- result$objfunc
}
# objective function by beta
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(betas,objfuncs,type="o",log="x")
table(result$block) # partition block of cv

(beta.best <- betas[which.min(objfuncs)])  
A <- nmfkc.kernel(U,beta=beta.best)
result <- nmfkc(Y,A,Q=1)
result$r.squared

# fitted curve
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(x,as.vector(Y))
lines(x,as.vector(result$XB),col=2,lwd=2)

# fitted curve for new data
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
V <- matrix(seq(from=min(U),to=max(U),length=100),ncol=100)
A <- nmfkc.kernel(U,V,beta=beta.best)
XB <- result$X %*% result$C %*% A
plot(x,as.vector(Y))
lines(V,as.vector(XB),col=2,lwd=2)
```

## 5. Growth curve model
- Orthodont
``` r
library(nmfkc)
library(nlme)
d <- Orthodont
head(d)
t <- unique(d$age)
Y <- matrix(d$distance,nrow=length(t))
colnames(Y) <- unique(d$Subject)
rownames(Y) <- t

Q <- 2
Male <- 1*(d$Sex=="Male")[d$age==8]
table(Male)
A <- rbind(rep(1,ncol(Y)),Male)
print(A)
result <- nmfkc(Y,A,Q=Q,epsilon=1e-8)
result$r.squared

# parameter matrix and coefficients by gender
result$C
(B <- t(unique(t(result$B))))

# individual fit
plot(t,Y[,1],ylim=range(Y),type="n")
mycol <- ifelse(Male==1,4,2)
for(n in 1:ncol(Y)){
  lines(t,Y[,n],col=mycol[n])
}
YHAT <- result$X %*% B
lines(t,YHAT[,1],col=4,lwd=5)
lines(t,YHAT[,2],col=2,lwd=5)
```

# Author
-  Kenichi Satoh, [homepage](https://sites.google.com/view/ksatoh/english)
