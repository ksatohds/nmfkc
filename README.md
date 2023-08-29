# Installation of **nmfkcreg** package

You can install the development version of nmfkcreg from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ksatohds/nmfkcreg")
```

# Functions

There are three functions in **nmfkcreg** package.

- **nmfkcreg** function is used for optimization of $X$ and $C$
- **nmfkcreg.cv** function is used for k-fold cross-validation
- **create.kernel** function is used for creating kernel matrix $A$ from covariate matrix $U$

# Statistical model

1. An ordinary LM (Linear Model) can be written as $y \approx Xb$, where $y$ is the observation, $X=(x_1,\cdots,x_Q)$ is the design matrix, and $b$ is a vector of regression coefficients. 

2. Using the common $X$ among $N$ individual observations, $[y_1,y_2,...,y_N] \approx X[b_1,...,b_N]$，i.e. $Y \approx XB$, this is called NMF (Non-Negative Matrix Factorisation). The difference between NMF and LM is that $X$ on NMF is optimized as well as $B$.

3. Since all the matrices are non-negative, the components of the regression coefficients $b_n$ for each $n$ are also non-negative. Therefore we know the contribution of the bases $x_1,.., x_Q$ and those percentages or probabilities $p_{1,n},...,p_{Q,n}$ can be used for clustering.

4. Furthermore, we consider that the regression coefficient matrix is described by the covariate matrix $A$, i.e. $B=CA$, where $C$ is the unknown parameter matrix. Here we propose to use a kernel matrix with high explanatory power as A.

# Matrices

The goal of **nmfkcreg** is to optimize $X$ and $C$ on the NMF (Non-negative Matrix Factorization) kernel covariate regression model $Y \approx X C A$ where $Y$ and $A$ are given.

- $Y(P,N)=(y_1,...y_N)$: **given** observation matrix
- $A(N,N)$: **given** kernel matrix of which $(i,j)$ element can be 
written as $K(u_i,u_j)=exp(−\beta|u_i-u_j|^2)$ here $U=(u_1,...u_N)$ 
is covariate matrix. Note that identity matrix is used as A=diag(ncol(Y)) in the case where there are no covariate. Or matrix $A(R,N)$ having N columns can be accepted.
- $X(P,Q)$: **unknown** basis matrix whose column sum is 1 and Q<=min{P,N}.
Q is the number of basis (rank).
-  $C(Q,N)$: **unkown** parameter matrix which is described by $\Theta$ in the paper Satoh (2023). 
- $B(Q,N)=CA$ is regression coefficient matrix. When $A$ is identity matrix, $B=C$. 
- $P(Q,N)$ is probability matrix for soft clustering based on regression coefficient matrix B. It is given by P <- t(t(B)/colSums(B)) whose column sum is 1.

# References

- Satoh, K. (2023) On Non-negative Matrix Factorization Using Gaussian Kernels as Covariates, Japanese Journal of Applied Statistics 52 (2), in press. (in Japanese)
[(preprint)](https://drive.google.com/file/d/1MnbJOPlcm0hn27WpP8rvcAgzy5X2E53B/view?usp=sharing)
- Satoh, K. (2022) Soft Clustering Based on Non-negative Matrix
Factorization for Longitudinal Data, Japanese Journal of Applied Statistics, 51 (1&2), 1-18. https://doi.org/10.5023/jappstat.51.1 (in Japanese)
- Ding, C., Tao, L., Wei, P. and Haesun, P. (2006)
Orthogonal Nonnegative Matrix Tri-Factorizations for Clustering,
Proceedings of the 12th ACM SIGKDD international conference on Knowledge discovery and data mining, 126-135.
- Potthoff, R. F. and Roy, S. N. (1964) A generalized multivariate analysis of variance model useful especially for growth curve problems,
Biometrika, 51, 313–326.


# References (in Japanese)
- 佐藤健一 (2023) ガウスカーネルを共変量に用いた非負値行列因子分解について, 応用統計学 52 (2), 印刷中. [(プレプリント)](https://drive.google.com/file/d/1MnbJOPlcm0hn27WpP8rvcAgzy5X2E53B/view?usp=sharing)
- 佐藤健一 (2022) 経時測定データに対する非負値行列因子分解によるソフトクラスタリングについて, 応用統計学, 51(1-2), 1-18. https://doi.org/10.5023/jappstat.51.1


# Examples

1. Dimension reduction: iris
2. Soft clulustering: basketball players and statistics
3. Hard clulustering: COVID-19 in Japan
4. Topic model: data_corpus_inaugural
5. Spatio-temporal Analysis: CanadianWeather
6. Comparison between NMF and ordinary LM: PimaIndiansDiabetes2
7. Kernel ridge regression: mcycle
8. Comparison between given covariate and kernel: Orthodont

## 1. Dimension reduction 
- iris
``` r
library(nmfkcreg)
Y <- t(iris[,-5])
result <- nmfkcreg(Y,Q=2)

# visualization of some results
plot(result$objfunc.iter) # convergence

# goodness of fit
plot(as.vector(result$XB),as.vector(Y),
main=paste0("r.squared=",round(result$r.squared,3)))
abline(a=0,b=1,col=2)

# dimension reduction based on regression coefficient B
labels <- as.numeric(iris[,5])
plot(t(result$B))
``` 

## 2. Soft clulustering
- basketball players and statistics
- https://rpubs.com/sirinya/847402
``` r
library(nmfkcreg)
d <- read.csv("http://datasets.flowingdata.com/ppg2008.csv")
y <- d[,-1]
rownames(y) <- d[,1]

Y <- t(y)
result <- nmfkcreg(Y,Q=2)

# visualization of some results
plot(result$objfunc.iter) # convergence

# goodness of fit
plot(as.vector(result$XB),as.vector(Y),
     main=paste0("r.squared=",round(result$r.squared,3)))
abline(a=0,b=1,col=2)

# individual fit
n <- 1
f <- rbind(Y[,n],result$XB[,n])
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

## 3. Hard clulustering
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
result <- nmfkcreg(Y,Q=7)
result$r.squared # goodness of fit

# basis function of which sum is 1
Q <- ncol(result$X)
par(mfrow=c(Q,1),mar=c(0,0,0,0))
for(q in 1:Q){
  barplot(result$X[,q],col=q+1,border=q+1,las=3,
    ylim=range(result$X),ylab=paste0("topic",q)) 
  legend("left",fill=q+1,legend=q)
}

# hard clulustering based on P
mycol <- result$cluster
library(NipponMap)
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
JapanPrefMap(col=mycol+1)
legend("left",fill=1:Q+1,legend=1:Q)
``` 

## 4. Topic model
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
U <- t(as.matrix(corp$Year))
Y <- t(d)
Q <- 3
result <- nmfkcreg(Y,Q=Q)
result$r.squared # coefficient of determination
colnames(result$P) <- corp$Year
barplot(result$P,col=1:Q+1,legend=T,las=3,ylab="Probability of topic")

#------------------
# with covariates using covariate matrix U
#------------------
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
plot(betas,objfuncs,type="o",log="x")
(best.beta <- betas[which.min(objfuncs)])

# create kernel with best beta
A <- create.kernel(U,beta=best.beta)
result <- nmfkcreg(Y,A,Q)
result$r.squared # less than nmf without covariates

# Topic probability changing over time
colnames(result$P) <- corp$Year
barplot(result$P,col=1:Q+1,legend=T,las=3,ylab="Probability of topic")

# frequent words that constitute each topic
par(mfrow=c(3,1))
for(q in 1:Q){
  f <- t(result$X[,q])
  names(f) <- rownames(result$X)
  index <- order(f,decreasing=T) 
  f <- f[index]
  barplot(f[1:30],col=q+1,las=3,ylab="probability",
    ylim=range(result$X),main=paste0("topic ",q))
}

# probability of topic at each frequent word
par(mfrow=c(1,1))
f <- prop.table(t(result$X[1:30,]),2)
barplot(f,col=1:Q+1,las=3,ylim=c(0,1),beside=T,legend=T,
  ylab="probability (relative contribution)")
abline(h=c(0.6,0.8),lty=3)
``` 

## 5. Spatio-temporal Analysis
-  CanadianWeather
``` r
library(nmfkcreg)
library(fda)
data(CanadianWeather)
d <- CanadianWeather$dailyAv[,,1]
Y <- d-min(d)
u0 <- CanadianWeather$coordinates[,2:1]
u0[,1] <- -u0[,1]
u <- t(u0)
umin <- apply(u,1,min)
umax <- apply(u,1,max)
U <- (u-umin)/(umax-umin) # normalization

#------------------
# without covariate
#------------------
result <- nmfkcreg(Y,Q=2)

# visualization of some results
plot(result$objfunc.iter) # convergence
result$r.squared # coefficient of determination

# individual fit
n <- 1
plot(Y[,n],main=colnames(Y)[n]) # observation
lines(result$XB[,n],col=2) # fitted values
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

#------------------
# with covariates using kernel matrix A
#------------------
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
plot(betas,objfuncs,type="o",log="x")
(best.beta <- betas[which.min(objfuncs)])

# create kernel with best beta
A <- create.kernel(U,beta=best.beta)
result <- nmfkcreg(Y,A,Q=2)
result$r.squared # less than nmf without covariates

# soft clulustering based on P by using covariates
library(akima)
q <- 2
result.interp <- interp(U[1,],U[2,],result$P[q,])
filled.contour(result.interp,
  xlab=rownames(U)[1],ylab=rownames(U)[2],
  plot.axes={
    points(t(U),col=3,pch=19)
    text(t(U),colnames(Y),pos=4)
  }
)
```

## 6. Comparison between NMF and ordinary LM
- PimaIndiansDiabetes2
``` r
library(nmfkcreg)
library(mlbench)
data(PimaIndiansDiabetes2)
d <- PimaIndiansDiabetes2
colnames(d)
d <- d[,-c(1,9)]
index <- complete.cases(d)
d <- d[index,]  # remove rows including NA's
res0 <- lm(glucose~.,d) # ordinary linear model

# preparation for NMF
# Y(1,N)~X(1,1)C(1,R)A(R,N) where X=1, Q=1 and R=7
Y <- t(as.matrix(d$glucose))
dim(Y) # 1*N
A <- t(as.matrix(cbind(1,d[,-1])))
dim(A) # R*N
library(nmfkcreg)
result <- nmfkcreg(Y,A,Q=1,epsilon=1e-15,maxit=20000)
plot(result$objfunc.iter,log="xy") # convergence
result$r.squared # coefficient of determination

# comparison between NMF and LM
f <- rbind(res0$coefficients,result$C)
rownames(f) <- c("LM","NMF")
print(f)
```

## 7. Kernel ridge regression
- mcycle
``` r
library(nmfkcreg)
library(MASS)
d <- mcycle
y0 <- d$accel
y <- y0-min(y0)
Y <- t(as.matrix(y))
U <- t(as.matrix(d$times))

# cv for optimization of beta and gamma
betas <- c(1,2,5,10,20)/100
gammas <- c(0,0.01,0.1)
bg <- expand.grid(betas,gammas)
objfuncs <- 0*(1:nrow(bg))
for(i in 1:nrow(bg)){
  print(i)
  A <- create.kernel(U,beta=bg[i,1])
  result <- nmfkcreg.cv(Y,A,gamma=bg[i,2],Q=1,div=10)
  objfuncs[i] <- result$objfunc
}
table(result$block) # partition block of cv
  
# objective function by beta and gamma
plot(1:nrow(bg),objfuncs,type="o")
text(1:nrow(bg),objfuncs,labels=paste0("b",bg[,1]),pos=3,col=2)
text(1:nrow(bg),objfuncs,labels=paste0("g",bg[,2]),pos=4,col=4)

(bg.best <- unlist(bg[which.min(objfuncs),]))  
A <- create.kernel(U,beta=bg.best[1])
result <- nmfkcreg(Y,A,Q=1,gamma=bg.best[2])

# visualization of some results
plot(result$objfunc.iter) # convergence  

# fitted curve
plot(d$times,as.vector(Y),
  main=paste0("r.squared=",round(result$r.squared,3)))
  lines(d$times,as.vector(result$XB),col=2)
```

## 8. Comparison between given covariate and kernel
- Orthodont
``` r
library(nlme)
d <- Orthodont
head(d)
Y <- matrix(d$distance,nrow=4)
N <- ncol(Y)
colnames(Y) <- unique(d$Subject)
t <- unique(d$age)
rownames(Y) <- t
Male8 <- 1*(d$Sex=="Male")[d$age==8]
U <- rbind(rep(1,N),Male8)

# ordinary covariate
library(nmfkcreg)
result.U <- nmfkcreg(Y,U,Q=2,epsilon=1e-8,maxit=20000)

# kernel covariate
betas <- c(0.05,0.1,0.2,0.5,1)
objfuncs <- 0*betas
for(i in 1:length(betas)){
  print(i)
  A <- create.kernel(U,beta=betas[i])
  result <- nmfkcreg(Y,A,Q=2,epsilon=1e-8,maxit=20000)
  objfuncs[i] <- result$objfunc
}
# objective function by beta
plot(betas,objfuncs,type="o",log="x")
(best.beta <- betas[which.min(objfuncs)])

A <- create.kernel(U,beta=best.beta)
result.A <- nmfkcreg(Y,A,Q=2,epsilon=1e-8,maxit=20000)
mycol <- ifelse(Male8==1,4,2)
matplot(t,Y,type="l",col=mycol,
        main=paste0("r.squared=",round(result.A$r.squared,5)))
matplot(t,result.U$XB,type="l",col=mycol,lwd=10,add=T)
matplot(t,result.A$XB,type="l",col=3,lwd=2,add=T)
mylab <- paste0("U(r.squared=",round(result.U$r.squared,5),")")
legend("topleft",legend=c("Male","Female",mylab),fill=c(4,2,3))
```

# Author
-  Kenichi SATOH, Ph.D.
-  Professor, [Faculty of Data Science](https://www.ds.shiga-u.ac.jp/), 
-  [Shiga University](https://www.shiga-u.ac.jp/), JAPAN
- [HomePage](https://sites.google.com/view/ksatoh/english)
