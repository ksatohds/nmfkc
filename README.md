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
- **nmfkc.rank** function is used to diagnose rank selection using the figure. Note that this function is under construction.

# Statistical model

1. An ordinary Linear Model (LM) can be written as $y \approx Xb$, where $y$ is the observation, $X=[x_1,...,x_Q]$ is the design matrix, and $b$ is a vector of regression coefficients. 

2. Using the common $X$ among $N$ individual observations, $[y_1,y_2,...,y_N] \approx X[b_1,...,b_N]$，i.e. $Y \approx XB$, this is called Non-Negative Matrix Factorization (NMF). The most significant difference between NMF and LM is that $X$ on NMF is optimized as well as $B$.

3. Since all the matrices are non-negative, the components of the regression coefficient $b=(b_1,...,b_Q)$ are also non-negative. Therefore the proportion $p_q=b_q/(b_1+...+b_Q)$, $q=1,...,Q$ can be used for soft clustering.

4. Furthermore, the coefficient $b$ can be explained by covariate $a$, i.e.,
$b=\Theta a$ or $B=\Theta A$ where $A=[a_1,a_2,...,a_N]$ and $\Theta$ is the parameter matrix to be optimized.
Satoh (2024) proposed Gaussian kernel function as covariates.
Formally the model is contained in tri-NMF by Ding et al. (2006) and the update formula for optimization has been derived. The model can be described as a product of three matrices and is related to the growth curve model by Potthoff and Roy (1964). 

# Matrices

The goal of **nmfkc** is to optimize $X(P,Q)$ and $C(Q,R)$ on the Non-negative Matrix Factorization, $Y(P,N) \approx X(P,Q)C(Q,R)A(R,N)=XB(Q,N)$ where $Y(P,N)$ and $A(R,N)$ are given.

- $Y(P,N)=(y_1,...y_N)$: **given** observation matrix
- $A(R,N)$: **given** covariate matrix.
  The covariate matrix can be created by nmfkc.kernel function.
  Note that identity matrix is used when there are no covariates.
- $X(P,Q)$: **unknown** basis matrix. Q is the number of basis (rank) and Q<=min(P,N).
-  $C(Q,R)$: **unknown** parameter matrix which is described by $\Theta$ in the paper Satoh (2024). 
- $B(Q,N)=C(Q,R)A(R,N)$ is coefficient matrix.

# Source

- Satoh, K. (2024) Applying Non-negative Matrix Factorization with Covariates to the Longitudinal Data as Growth Curve Model. arXiv preprint arXiv:2403.05359. https://arxiv.org/abs/2403.05359

# References

- Ding, C., Li, T., Peng, W. and Park, H. (2006) Orthogonal Nonnegative Matrix Tri-Factorizations for Clustering, Proceedings of the 12th ACM SIGKDD international conference on Knowledge discovery and data mining, 126-135. https://doi.org/10.1145/1150402.1150420
- Potthoff, R.F., and Roy, S.N. (1964). A generalized multivariate analysis of variance model useful especially for growth curve problems. Biometrika, 51, 313-326. https://doi.org/10.2307/2334137

# Examples

1. Longitudinal data: COVID-19 in Japan
2. Spatiotemporal Analysis: CanadianWeather
3. Topic model: data_corpus_inaugural
4. Origin-Destination (OD) data: Japanese Inter-prefecture flow
5. Kernel ridge regression: mcycle
6. Growth curve model: Orthodont
7. Binary repeated measures: Table 6, Koch et al.(1977)

## 1. Longitudinal data
- COVID-19 in Japan
- https://www3.nhk.or.jp/news/special/coronavirus/data/
``` r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("NipponMap")

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

# rank selection diagnostics
library(nmfkc)
nmfkc.rank(Y,Q=2:8)

# nmf
Q <- 4
result <- nmfkc(Y,Q=Q)
result$r.squared # goodness of fit

# individual fit
par(mfrow=c(7,7),mar=c(0,0,0,0)+0.1,cex=1)
for(n in 1:ncol(Y)){
  plot(Y[,n],axes=F,type="l") # observation
  lines(result$XB[,n],col=2) # fitted values
  legend("topleft",legend=colnames(Y)[n],x.intersp=-0.5,bty="n")  
  box()
}

# basis function of which sum is 1
par(mfrow=c(Q,1),mar=c(0,0,0,0),cex=1)
for(q in 1:Q){
  barplot(result$X[,q],col=q+1,border=q+1,las=3,
    ylim=range(result$X),ylab=paste0("topic ",q)) 
  legend("left",fill=q+1,legend=q)
}

# cluster membership probability based on coefficients
n <- 1
result$B[,n]
result$B[,n]/sum(result$B[,n])
result$B.prob[,n]

# soft clustering based on B.prob
library(NipponMap)
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
jmap <- JapanPrefMap(col="white",axes=TRUE)
stars(x=t(result$B.prob),scale=F,
      locations=jmap,key.loc =c(145,34),
      draw.segments=T,len=0.7,labels=NULL,
      col.segments=c(1:Q)+1,add=T)

# heatmap
heatmap(t(result$B.prob))

# hard clustering based on B.cluster
table(result$B.cluster)
``` 

## 2. Spatiotemporal Analysis
-  CanadianWeather
``` r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("fda")

library(fda)
data(CanadianWeather)
d <- CanadianWeather$dailyAv[,,1]
Y <- d-min(d)
u0 <- CanadianWeather$coordinates[,2:1]
u0[,1] <- -u0[,1]

#------------------
# without covariate
#------------------
library(nmfkc)
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
  legend("topleft",legend=colnames(Y)[n],x.intersp=-0.5,bty="n")  
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
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("quanteda")

#------------------
# text analysis
#------------------
library(quanteda)
corp <- corpus(data_corpus_inaugural)
head(corp,3)
tail(corp,3)
tok <- tokens(corp)
tok <- tokens_remove(tok,pattern=stopwords("en",source="snowball"))
df <- dfm(tok)
df <- dfm_select(df,min_nchar=3)
df <- dfm_trim(df,min_termfreq=100)
d <- as.matrix(df)
index <- order(colSums(d),decreasing=T) 
d <- d[,index] # document-word matrix
paste0(colnames(d)[1:30],"(",colSums(d)[1:30],")") # Top 30 most frequent words

#------------------
# without covariates
#------------------
Y <- t(d)
Y[1:20,c(1,ncol(Y))]
Q <- 3
library(nmfkc)
result <- nmfkc(Y,Q=Q)
result$r.squared # coefficient of determination

# topic probabilities
par(mfrow=c(1,1),mar=c(10,4,4,2)+0.1,cex=1)
barplot(result$B.prob,col=1:Q+1,legend=T,las=3,
  ylab="Probabilities of topics")

# basis function of which sum is 1
par(mfrow=c(Q,1),mar=c(8,4,1,1),cex=0.6)
for(q in 1:Q){
  barplot(result$X[,q],col=q+1,border=q+1,las=3,
    ylim=range(result$X),ylab=paste0("topic ",q)) 
}

# contribution of words to each topics
Xp <- prop.table(result$X,1)
par(mfrow=c(1,1),mar=c(8,4,1,1)+0.1,cex=0.6)
barplot(t(Xp),las=3,col=1:Q+1,
  ylab="Proportion of words on each topic")
legend("topleft",fill=1:Q+1,legend=paste0("topic ",1:Q))
for(q in 1:Q){
  message(paste0("----- featured words on topic [",q,"] -----"))
  print(paste0(rownames(Xp),"(",rowSums(Y),")",round(100*Xp[,q],1),"%")[Xp[,q]>=0.5])
}

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

## 4. Origin-Destination (OD) data: Japanese Inter-prefecture flow
- e-stat run by Japanese Government Statistics
``` r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("httr")
# install.packages("readxl")
# install.packages("alluvial")
# install.packages("NipponMap")
# install.packages("RColorBrewer")

# download data and its formatting
library(httr)
# https://www.e-stat.go.jp/stat-search/files?stat_infid=000040170612
url <- "https://www.e-stat.go.jp/stat-search/file-download?statInfId=000040170612&fileKind=0"
GET(url, write_disk(tmp <- tempfile(fileext=".xlsx")))
library(readxl)
d <- as.data.frame(read_xlsx(tmp,sheet=1,skip=2))
colnames(d)
d <- d[d[,1]=="1",]
d <- d[d[,3]=="02",]
d <- d[d[,5]!="100",]
d <- d[d[,7]!="100",]
pref <- unique(d[,6])

#------------------
# without covariates
#------------------
Y <- matrix(NA,nrow=47,ncol=47)
colnames(Y) <- pref
rownames(Y) <- pref
d[,5] <- as.numeric(d[,5])
d[,7] <- as.numeric(d[,7])
for(i in 1:47)for(j in 1:47){
  Y[i,j] <- d[which(d[,5]==i&d[,7]==j),9]
}
Y <- log(1+Y)
Y[1:6,1:6]

# rank selection diagnostics
library(nmfkc)
nmfkc.rank(Y,Q=2:12)

# nmf
Q0 <- 7
res <- nmfkc(Y,Q=Q0)
plot(res$objfunc.iter,type="o",
     main=paste0("Q=",Q0,", R^2=",round(res$r.squared,3)))
     
# Silhouette: intra-cluster accumulation and inter-cluster discrepancy
si <- res$criterion$silhouette
barplot(si$silhouette,horiz=T,las=1,col=si$cluster+1,cex.names=0.5,xlab="")
abline(v=si$silhouette.mean,lty=3)
legend("bottomleft",fill=1:Q0,legend=1:Q0)

# Sankey diagram: stability of hard clustering by 
Q <- 6:8
cluster <- NULL
for(i in 1:length(Q)){
  res <- nmfkc(Y,Q=Q[i])
  cluster <- cbind(cluster,res$B.cluster)
}
library(alluvial)
alluvial(cluster,freq=1,axis_labels=paste0("Q=",Q),cex=2,
         col=cluster[,2]+1,border=cluster[,2]+1)

# basis function of which sum is 1
library(NipponMap)
library(RColorBrewer)
mypalette <- brewer.pal(9,"YlOrRd")
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1,cex=0.6)
for(j in 1:Q0){
  cutp <- as.numeric(
    cut(res$B.prob[j,],
        breaks=seq(from=0,to=1,length=10),
        include.lowest=T))
  mycol <- mypalette[cutp]
  par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
  jmap <- JapanPrefMap(col=mycol,axes=TRUE,
                       main=paste0("basis[",j,"]"))
  text(jmap,pref,cex=0.5)  
}

# soft clustering based on B.prob
mypalette <- brewer.pal(12,"Paired")
tp <- t(res$B.prob)
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
jmap <- JapanPrefMap(col="white",axes=TRUE)
stars(x=tp,scale=F,locations=jmap,key.loc =c(145,34),
      draw.segments=TRUE,len=1,labels=NULL,
      col.segments=mypalette[1:Q0],add=T)
title(main="Inter-prefecture flow: weekdays - operations")

# hard clustering based on B.cluster
par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
jmap <- JapanPrefMap(col=mypalette[res$B.cluster],axes=TRUE)
text(jmap,pref,cex=0.5)  
legend("topleft",fill=mypalette[1:Q0],
       legend=1:Q0,title="basis")
title(main="Inter-prefecture flow: weekdays - operations")
```

## 5. Kernel ridge regression
- mcycle
``` r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")

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
library(nmfkc)
result <- nmfkc(Y,A,Q=1)
result$r.squared
par(mfrow=c(1,1),mar=c(5,4,2,2)+0.1,cex=1)
plot(U,Y)
lines(U,result$XB,col=2)

# cv for optimization of beta
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

# Periodic kernel function
A <- nmfkc.kernel(U,method="Periodic",beta=c(11,0.07))
result <- nmfkc(Y,A,Q=1)
A <- nmfkc.kernel(U,V,method="Periodic",beta=c(11,0.07))
XB <- result$X %*% result$C %*% A
lines(as.vector(V),as.vector(XB),col=4,lwd=2)
legend("topleft",legend=c("Gaussian","Periodic"),fill=c(2,4))
```

## 6. Growth curve model
- Orthodont
``` r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# install.packages("nlme")

library(nlme)
d <- Orthodont
head(d)
t <- unique(d$age)
Y <- matrix(d$distance,nrow=length(t))
colnames(Y) <- unique(d$Subject)
rownames(Y) <- t

#------------------
# with covariates
#------------------
Male <- 1*(d$Sex=="Male")[d$age==8]
table(Male)
A <- rbind(rep(1,ncol(Y)),Male)
rownames(A) <- c("Const","Male")

# nmf with covariates
library(nmfkc)
result <- nmfkc(Y,A,Q=2,epsilon=1e-8)
result$r.squared

# basis matrix
print.table(round(result$X,2),zero.print="")

# unique covariates and coefficient matrix
(A0 <- t(unique(t(A))))
B <- result$C %*% A0

# individual fit
plot(t,Y[,1],ylim=range(Y),type="n")
mycol <- ifelse(Male==1,4,2)
for(n in 1:ncol(Y)){
  lines(t,Y[,n],col=mycol[n])
}
XB <- result$X %*% B
lines(t,XB[,1],col=4,lwd=5)
lines(t,XB[,2],col=2,lwd=5)
legend("topleft",title="Sex",legend=c("Male","Female"),fill=c(4,2))
```

## 7. Binary repeated measures
- Table 6, Koch et al.(1977) Biometrics, 33(1), 133–158. https://doi.org/10.2307/2529309
- Repeated categorical outcome analysis, https://wnarifin.github.io/medstat.html
``` r
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")

Id <- rep(1:340,each=3)
Mild <- rep(c(1,0),c(150*3,190*3))
NewDrug <- rep(c(0,1,0,1),c(80*3,70*3,100*3,90*3))
Week <- rep(c(1,2,4),times=340)
Normal <- unlist(c(
  rep(list(c(1,1,1),c(1,1,0),c(1,0,1),c(1,0,0),c(0,1,1),c(0,1,0),c(0,0,1),c(0,0,0)),c(16,13,9,3,14,4,15,6)),
  rep(list(c(1,1,1),c(1,1,0),c(1,0,1),c(1,0,0),c(0,1,1),c(0,1,0),c(0,0,1),c(0,0,0)),c(31,0,6,0,22,2,9,0)),
  rep(list(c(1,1,1),c(1,1,0),c(1,0,1),c(1,0,0),c(0,1,1),c(0,1,0),c(0,0,1),c(0,0,0)),c(2,2,8,9,9,15,27,28)),
  rep(list(c(1,1,1),c(1,1,0),c(1,0,1),c(1,0,0),c(0,1,1),c(0,1,0),c(0,0,1),c(0,0,0)),c(7,2,5,2,31,5,32,6))
))
d <- data.frame(Id,Mild,NewDrug,Week,Normal)
ftable(xtabs(~Mild+NewDrug+Week+Normal))

# observation matrix and covariate matrix
Y <- matrix(d$Normal,nrow=3,ncol=340)
A <- matrix(0,nrow=3,ncol=340)
A[1,] <- 1
A[2,] <- d$Mild[d$Week==1]
A[3,] <- d$NewDrug[d$Week==1]
rownames(A) <- c("Const","Mild","NewDrug")

# nmf with covariates
library(nmfkc)
result <- nmfkc(Y,A,Q=2,epsilon=1e-8)
plot(result)

# unique covariates and coefficient matrix
(A0 <- t(unique(t(A))))
B <- result$C %*% A0

# fitted values for coefficient matrix
t <- unique(d$Week)
plot(t,Y[,1],ylim=range(Y),type="n",ylab="Normal",xlab="Week")
XB <- result$X %*% B
mycol <- c(2,7,4,5)
for(j in 1:ncol(XB))lines(t,XB[,j],col=mycol[j],lwd=5)
legend("bottomright",
       title="Diagnosis & Treatment",
       legend=c("Mild & Standard","Mild & NewDrug","Sever & Standard","Sever & NewDrug"),
       fill=mycol)
```


# Author
-  Kenichi Satoh, [homepage](https://sites.google.com/view/ksatoh/english)
