.onAttach <- function(...) {
  packageStartupMessage("Last update on 15 FEB 2025")
  packageStartupMessage("https://github.com/ksatohds/nmfkc")
}

#' @title Creating observation and covariates for the vector autoregressive model
#' @description \code{nmfkc.ar} create observation matrix and covariate matrix according to the degree of the autoregressive model
#' @param Y Observation matrices with columns in ascending order of measurement time point
#' @param degree lag degree for the autoregressive model, and the default value is 1.
#' @param intercept The default value is TRUE to add the intercept to the covariate matrix
#' @return Y: observation matrix according to the degree of the autoregressive model
#' @return A: covariate matrix according to the degree of the autoregressive model
#' @return A.columns: subscript matrix used to create A
#' @return degree.max: 10log10(N) according to ar function in stats package
#' @export
#' @source Satoh, K.(2025) Applying non-negative Matrix Factorization with Covariates to Multivariate Time Series Data as a Vector Autoregression Model. arXiv:2501.17446. \url{https://arxiv.org/abs/2501.17446}
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' d <- AirPassengers
#' time <- time(ts(1:length(d),start=c(1949,1),frequency=12))
#' time.vec <- round(as.vector(t(time)),2)
#' Y0 <- matrix(as.vector(d),nrow=1)
#' colnames(Y0) <- time.vec
#' rownames(Y0) <- "t"
#' # nmf with covariates
#' Q <- 1
#' D <- 12
#' a <- nmfkc.ar(Y0,degree=D,intercept=TRUE); Y <- a$Y; A <- a$A
#' res <- nmfkc(Y=Y,A=A,Q=Q,prefix="Factor",epsilon=1e-6)
#' res$r.squared
#' # coefficients
#' print.table(round(res$C,2),zero.print="")
#' # fitted curve
#' plot(as.numeric(colnames(Y)),as.vector(Y),type="l",col=1,xlab="",ylab="AirPassengers")
#' lines(as.numeric(colnames(Y)),as.vector(res$XB),col=2)

nmfkc.ar <- function(Y,degree=1,intercept=T){
  if(is.vector(Y)) Y <- matrix(Y,nrow=1)
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  N <- ncol(Y)
  A.columns <- NULL
  for(i in degree:1)A.columns <- rbind(A.columns,i:(i+ncol(Y)-degree-1))
  A <- NULL
  for(i in 1:nrow(A.columns))A <- rbind(A,Y[,A.columns[i,]])
  if(is.null(rownames(Y))) rownames(Y) <- 1:nrow(Y)
  label <- NULL
  for(i in 1:degree)label <- c(label,paste0(rownames(Y),"_",i))
  if(intercept){
    A <- rbind(A,1)
    rownames(A) <- c(label,"(Intercept)")
  }
  Ya <- Y[,A.columns[1,]+1]
  if(!is.matrix(Ya)){
    Ya <- matrix(Ya,nrow=1)
    colnames(Ya) <- colnames(Y)[A.columns[1,]+1]
    rownames(Ya) <- rownames(Y)[1]
  }
  degree.max <- min(ncol(Ya),floor(10*log10(ncol(Ya))))
  list(Y=Ya,A=A,A.columns=A.columns)
}


#' @title DOT language script for the vector autoregressive model
#' @description \code{nmfkc.ar.DOT} create scripts in DOT language function
#' @param x return value of nmfkc function for vector autoregressive
#' @param degree max lag degree to visualize, and the default value is 1.
#' @param intercept The default value is FALSE to add
#' @param digits integer indicating the number of decimal places
#' @param threshold display parameters that exceed the threshold value
#' @param rankdir Direction of node placement based on DOT language, default is RL, LR, TB, BT can also be specified.
#' @return scripts for dot function of DOT package
#' @export

nmfkc.ar.DOT <- function(x,degree=1,intercept=FALSE,digits=1,threshold=10^(-digits),rankdir="RL"){
  X <- x$X; C <- round(x$C,digits); D <- min(ncol(C),degree)
  rownames(X) <- gsub(".","",rownames(X),fixed=T)
  colnames(X) <- gsub(".","",colnames(X),fixed=T)
  rownames(C) <- gsub(".","",rownames(C),fixed=T)
  colnames(C) <- gsub(".","",colnames(C),fixed=T)
  # rankdir=RL # rankdir=TB
  scr <- paste0('digraph XCA {graph [rankdir=',rankdir,' compound=true]; \n')
  # Y
  st <- 'subgraph cluster_Y{label="T" style="rounded"; \n'
  for(j in 1:nrow(X))st <- paste0(st,sprintf('  %s [shape=box]; \n',rownames(X)[j]))
  st <- paste0(st,'}; \n'); scr <- paste0(scr,st)
  # X and element
  st <- 'subgraph cluster_X{label="Latant Variables" style="rounded"; \n'
  for(j in 1:ncol(X))st <- paste0(st,sprintf('  %s [shape=ellipse]; \n',colnames(X)[j]))
  st <- paste0(st,'}; \n'); scr <- paste0(scr,st)
  # edge: X to Y
  for(i in 1:nrow(X))for(j in 1:ncol(X)){
    if(X[i,j]>=threshold){
      st <- sprintf(paste0('%s -> %s [label="%.',digits,'f"]; \n'),
              colnames(X)[j],rownames(X)[i],X[i,j]); scr <- paste0(scr,st)}
  }
  # intercept
  if(intercept==TRUE){
    for(i in 1:ncol(X))
      if(C[i,ncol(C)]>=threshold){
        st <- sprintf(paste0('Const%d [shape=circle label="%.',digits,'f"]; '),i,C[i,ncol(C)])
        st <- paste0(st,sprintf(paste0('Const%d -> %s; \n'),i,colnames(X)[i]))
        scr <- paste0(scr,st)
      }
  }
  # edge: T-k to X
  klist <- NULL; ktoplist <- NULL
  for(k in 1:D){
    Ck <- C[,(k-1)*nrow(X)+1:nrow(X)]
    if(is.matrix(Ck)==FALSE){
      Ck <- matrix(Ck,nrow=1)
      colnames(Ck) <- colnames(C)[k]
      rownames(Ck) <- colnames(X)[1]
    }
    if(max(Ck)>=threshold){
      klist <- c(klist,k)
      st <- sprintf('subgraph cluster_C%d{label="T-%d" style="rounded"; \n',k,k)
      ktop <- NULL
      for(j in 1:ncol(Ck)){
        if(max(Ck[,j])>=threshold){
          if(is.null(ktop)==TRUE){
            ktop <- j
            ktoplist <- c(ktoplist,ktop)
          }
          st <- paste0(st,sprintf('  %s [label="%s",shape=box]; \n',
                                  colnames(Ck)[j],rownames(X)[j]))
        }
      }
      st <- paste0(st,"}; \n");scr <- paste0(scr,st)
    }
    for(i in 1:nrow(Ck))for(j in 1:ncol(Ck)){
      if(Ck[i,j]>=threshold){
        st <- sprintf(paste0('%s -> %s [label="%.',digits,'f"]; \n'),
                      colnames(Ck)[j],rownames(Ck)[i],Ck[i,j])
        scr <- paste0(scr,st)
      }
    }
  }
  if(length(klist)>=2){
    for(k in 2:length(klist)){
      st <- sprintf('%s_%d -> %s_%d [ltail=cluster_C%d lhead=cluster_C%d style=invis]; \n',
              rownames(X)[ktoplist[k]],klist[k],
              rownames(X)[ktoplist[k-1]],klist[k-1],
                          klist[k],klist[k-1])
      scr <- paste0(scr,st)
    }
  }
  scr <- paste0(scr,"} \n")
  return(scr)
}


#' @title Optimizing degree for the autoregressive model
#' @description \code{nmfkc.ar.degree.cv} apply cross validation method for degree
#' @param Y observation matrix
#' @param Q rank of basis matrix and Q<=min(P,N) where Y(P,N)
#' @param degree degree vector for cv
#' @param intercept The default value is TRUE to add the intercept to the covariate matrix
#' @param div number of partition usually described as "k" of k-fold
#' @param seed integer used as argument in set.seed function
#' @param plot The default is plot=TRUE and draw a graph.
#' @return degree: best degree minimizes objective function
#' @return degree.max: maximum recommended degree in ar model
#' @return objfunc: objective functions
#' @export
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' d <- AirPassengers
#' time <- time(ts(1:length(d),start=c(1949,1),frequency=12))
#' time.vec <- round(as.vector(t(time)),2)
#' Y0 <- matrix(as.vector(d),nrow=1)
#' colnames(Y0) <- time.vec
#' rownames(Y0) <- "t"
#' # selection of degree
#' nmfkc.ar.degree.cv(Y=Y0,Q=1,degree=11:14)

nmfkc.ar.degree.cv <- function(Y,Q=2,degree=1:2,intercept=T,div=5,seed=123,plot=TRUE){
  objfuncs <- 0*(1:length(degree))
  for(i in 1:length(degree)){
    start.time <- Sys.time()
    packageStartupMessage(paste0("degree=",degree[i],"..."),appendLF=FALSE)
    a <- nmfkc.ar(Y=Y,degree=degree[i],intercept=intercept)
    result.cv <- nmfkc.cv(Y=a$Y,A=a$A,Q=Q,div=div,seed=seed)
    objfuncs[i] <- result.cv$objfunc/ncol(a$Y)
    end.time <- Sys.time()
    diff.time <- difftime(end.time,start.time,units="sec")
    diff.time.st <- ifelse(diff.time<=180,paste0(round(diff.time,1),"sec"),
                           paste0(round(diff.time/60,1),"min"))
    packageStartupMessage(diff.time.st)
  }
  i0 <- which.min(objfuncs)
  best.degree <- degree[i0]
  degree.max <- min(ncol(Y),floor(10*log10(ncol(Y))))
  if(plot){
    plot(degree,objfuncs,type="l",col=2,xlab=paste0("degree (max=",degree.max,")"),ylab="objfunc")
    graphics::points(degree[i0],objfuncs[i0],cex=3,col=2)
    graphics::text(degree,objfuncs,degree)
  }
  names(objfuncs) <- degree
  result <- list(degree=best.degree,degree.max=degree.max,objfunc=objfuncs)
  return(result)
}


#' @title Creating kernel matrix from covariates
#' @description \code{nmfkc.kernel} create kernel matrix from covariate matrix
#' @param U covariate matrix U(K,N)=(u_1,...,u_N) each row might be normalized in advance
#' @param V covariate matrix V(K,M)=(v_1,...,v_M) usually used for prediction, and if it is NULL, the default value is U.
#' @param beta The default parameter of kernel function is 0.5.
#' @param kernel The default kernel function is Gaussian kernel. For other functions, check by typing "nmfkc.kernel".
#' @param degree The default parameter of kernel function is 2.
#' @return kernel matrix A(N,M)
#' @export
#' @source Satoh, K. (2024) Applying Non-negative Matrix Factorization with Covariates to the Longitudinal Data as Growth Curve Model. arXiv preprint arXiv:2403.05359. \url{https://arxiv.org/abs/2403.05359}
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' Y <- matrix(cars$dist,nrow=1)
#' U <- matrix(c(5,10,15,20,25),nrow=1)
#' V <- matrix(cars$speed,nrow=1)
#' A <- nmfkc.kernel(U,V,beta=28/1000)
#' dim(A)
#' result <- nmfkc(Y,A,Q=1)
#' plot(as.vector(V),as.vector(Y))
#' lines(as.vector(V),as.vector(result$XB),col=2,lwd=2)

nmfkc.kernel <- function(U,V=NULL,beta=0.5,kernel="Gaussian",degree=2){
  if(is.null(V)==TRUE) V <- U
  kvec <- function(m){
    vm <- t(rep(1,ncol(U)) %o% V[,m])
    d <- colSums((U-vm)^2)^0.5
    k <- 0
    if(kernel=="Gaussian") k <- exp(-beta*d^2) # Gaussian
    if(kernel=="Exponential") k <- exp(-beta*d)
    if(kernel=="Periodic") k <- exp(-beta[1]*sin(beta[2]*d)^2)
    if(kernel=="Linear") k <- t(U) %*% V[,m]
    if(kernel=="NormalizedLinear") k <- diag(1/colSums(U^2)^0.5) %*% t(U) %*% V[,m]/sum(V[,m]^2)^0.5
    if(kernel=="Polynomial") k <- (t(U) %*% V[,m]+beta)^degree
    return(k)}
  A <- NULL; for(m in 1:ncol(V)) A <- cbind(A,kvec(m))
  if(min(A)<0){
    warning("The constructed matrix is not non-negative.")
  }
  return(A)
}

#' @title Optimizing beta of Gauss kernel function by cross validation method
#' @description \code{nmfkc.kernel.beta.cv} apply cross validation method for beta
#' @param Y observation matrix
#' @param Q rank of basis matrix and Q<=min(P,N) where Y(P,N)
#' @param U covariate matrix U(K,N)=(u_1,...,u_N) each row might be normalized in advance
#' @param V covariate matrix V(K,M)=(v_1,...,v_M) usually used for prediction, and if it is NULL, the default value is U.
#' @param beta parameter vector of kernel function for cv
#' @param kernel The default kernel function is Gaussian kernel. For other functions, check by typing "nmfkc.kernel".
#' @param degree The default parameter of kernel function is 2.
#' @param div number of partition usually described as "k" of k-fold
#' @param seed integer used as argument in set.seed function
#' @param plot The default is plot=TRUE and draw a graph.
#' @return beta: best parameter minimizes objective function
#' @return objfunc: objective functions
#' @export
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' Y <- matrix(cars$dist,nrow=1)
#' U <- matrix(c(5,10,15,20,25),nrow=1)
#' V <- matrix(cars$speed,nrow=1)
#' nmfkc.kernel.beta.cv(Y,Q=1,U,V,beta=25:30/1000)
#' A <- nmfkc.kernel(U,V,beta=28/1000)
#' result <- nmfkc(Y,A,Q=1)
#' plot(as.vector(V),as.vector(Y))
#' lines(as.vector(V),as.vector(result$XB),col=2,lwd=2)

nmfkc.kernel.beta.cv <- function(Y,Q=2,U,V=NULL,beta=c(0.1,0.2,0.5,1,2,5,10,20,50),
                                 kernel="Gaussian",degree=2,div=5,seed=123,plot=TRUE){
  objfuncs <- 0*(1:length(beta))
  for(i in 1:length(beta)){
    start.time <- Sys.time()
    packageStartupMessage(paste0("beta=",beta[i],"..."),appendLF=FALSE)
    A <- nmfkc.kernel(U=U,V=V,beta=beta[i],kernel=kernel,degree=degree)
    result <- nmfkc.cv(Y=Y,A=A,Q=Q,div=div,seed=seed)
    objfuncs[i] <- result$objfunc
    end.time <- Sys.time()
    diff.time <- difftime(end.time,start.time,units="sec")
    diff.time.st <- ifelse(diff.time<=180,paste0(round(diff.time,1),"sec"),
                           paste0(round(diff.time/60,1),"min"))
    packageStartupMessage(diff.time.st)
  }
  i0 <- which.min(objfuncs)
  beta.best <- beta[i0]
  if(plot){
    plot(beta,objfuncs,type="l",col=2,xlab="beta",ylab="objfunc",log="x")
    graphics::points(beta[i0],objfuncs[i0],cex=3,col=2)
    graphics::text(beta,objfuncs,beta)
  }
  names(objfuncs) <- beta
  result <- list(beta=beta.best,objfunc=objfuncs)
  return(result)
}


#' @title Optimizing NMF (Non-negative Matrix Factorization) with Kernel Covariate
#' @description \code{nmfkc} The goal of the package is to perform NMF (Non-negative Matrix Factorization) with Kernel Covariate described by Y~XCA=XB
#'  where observation matrix Y(P,N),
#'  covariate matrix A(R,N),
#'  basis matrix X(P,Q) and Q<=min(P,N),
#'  parameter matrix C(Q,R)
#'  and coefficient matrix B(Q,N).
#'  Note that Y(N,P) and A(R,N) are given, and X(P,Q) and C(Q,R) are unknown.
#' @param Y observation matrix
#' @param A covariate matrix. The default is NULL if without covariate.
#' @param Q rank of basis matrix and Q<=min(P,N)
#' @param gamma penalty parameter for parameter matrix C
#' @param epsilon positive convergence tolerance
#' @param maxit maximum number of iterations
#' @param method The default objective function is Euclid distance "EU", otherwise Kullback–Leibler divergence "KL"
#' @param X.restriction The default is X.restriction="colSums" and the column sum of basis matrix is 1, and it is interpreted as probability. The column of basis matrix is unit vector when X.restriction="colSqSums".
#' @param nstart The default is one. It is the "nstart" option of "kmeans" function used for the initial values of basis matrix.
#' @param seed integer used as argument in set.seed function
#' @param prefix prefix label for column labels of basis matrix and row labels of coefficient matrix
#' @param print.trace display current iteration every 10 times if print.trace=TRUE
#' @param print.dims display dimensions of matrix sizes if print.dim=TRUE. The default is set by  print.dim=FALSE.
#' @param save.time The default is TRUE. Some return values including CPCC and silhouette are skipped to save the computation time.
#' @param save.memory The default is FALSE. If save.memory=TRUE, save.time=TRUE and only the minimum necessary calculations are performed to save memory.
#' @return X: basis matrix. The column sum depends on X.restriction.
#' @return B: coefficient matrix, B=CA
#' @return XB: fitted values for observation matrix Y
#' @return C: parameter matrix
#' @return B.prob: probability matrix for soft clustering based on column vector of coefficient matrix B.
#' @return B.cluster: the number of the basis that takes the maximum value of each column of B.prob for hard clustering
#' @return X.prob: probability matrix for soft clustering based on row vector of basis matrix X.
#' @return X.cluster: the number of the basis that takes the maximum value of each row of X.prob for hard clustering
#' @return objfunc: last objective function
#' @return objfunc.iter: objective function at each iteration
#' @return r.squared: coefficient of determination R^2, squared correlation between Y and XB
#' @return criterion: several criteria for selecting rank including ICp, CPCC, silhouette
#' @export
#' @source Satoh, K. (2024) Applying Non-negative Matrix Factorization with Covariates to the Longitudinal Data as Growth Curve Model. arXiv preprint arXiv:2403.05359. \url{https://arxiv.org/abs/2403.05359}
#' @references Ding, C., Li, T., Peng, W. and Park, H. (2006) Orthogonal Nonnegative Matrix Tri-Factorizations for Clustering, Proceedings of the 12th ACM SIGKDD international conference on Knowledge discovery and data mining, 126-135. \url{https://doi.org/10.1145/1150402.1150420}
#' @references Potthoff, R.F. and Roy, S.N. (1964). A generalized multivariate analysis of variance model useful especially for growth curve problems. Biometrika, 51, 313-326. \url{https://doi.org/10.2307/2334137}
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example 1.
#' library(nmfkc)
#' X <- cbind(c(1,0,1),c(0,1,0))
#' B <- cbind(c(1,0),c(0,1),c(1,1))
#' Y <- X %*% B
#' rownames(Y) <- paste0("P",1:nrow(Y))
#' colnames(Y) <- paste0("N",1:ncol(Y))
#' print(X); print(B); print(Y)
#' library(nmfkc)
#' res <- nmfkc(Y,Q=2)
#' res$X
#' res$B
#'
#' # Example 2.
#' Y <- matrix(cars$dist,nrow=1)
#' A <- rbind(1,cars$speed)
#' result <- nmfkc(Y,A,Q=1)
#' plot(cars$speed,as.vector(Y))
#' lines(cars$speed,as.vector(result$XB),col=2,lwd=2)

nmfkc <- function(Y,A=NULL,Q=2,gamma=0,epsilon=1e-4,maxit=5000,method="EU",
                  X.restriction="colSums",nstart=1,seed=123,
                  prefix="Basis",print.trace=FALSE,print.dims=TRUE,save.time=TRUE,save.memory=FALSE){
  z <- function(x){
    x[is.nan(x)] <- 0
    x[is.infinite(x)] <- 0
    return(x)
  }
  mysilhouette <- function(B.prob,B.cluster){
    if(is.matrix(B.prob)){Q <- nrow(B.prob)}else{Q <- 1}
    if(Q==1){
      return(list(cluster=NA,silhouette=NA,
                  silhouette.means=NA,silhouette.mean=NA))
    }else{
      index <-!is.na(B.cluster)
      B.prob <- B.prob[,index]
      B.cluster <- B.cluster[index]
      cluster.means <- matrix(0,nrow=Q,ncol=Q)
      ns <- NULL
      cluster.list <- NULL
      for(q in 1:Q){
        ns <- c(ns,sum(B.cluster==q))
        cluster.list <- c(cluster.list,list(which(B.cluster==q)))
        cluster.means[,q] <- rowMeans(B.prob[,B.cluster==q,drop=F])
      }
      si <- 0*B.cluster
      neighbor.cluster <- 0*B.cluster
      for(q in 1:Q){
        for(i in cluster.list[[q]]){
          di <- colSums((cluster.means-B.prob[,i])^2)
          qn <- ifelse(order(di)[1]==q,order(di)[2],order(di)[1])
          neighbor.cluster[i] <- qn
          if(ns[q]==1){
            si[i] <- 0
          }else{
            ai <- sum(colSums((B.prob[,cluster.list[[q]],drop=F]-B.prob[,i])^2)^0.5)/(ns[q]-1)
            bi <- sum(colSums((B.prob[,cluster.list[[qn]],drop=F]-B.prob[,i])^2)^0.5)/ns[qn]
            si[i] <- (bi-ai)/max(ai,bi)
          }
        }
      }
      si.mean <- mean(si)
      si.sort.cluster.means <- 0*ns
      for(q in 1:Q){
        si.sort.cluster.means[q] <- mean(cluster.list[[q]])
      }
      si.sort <- NULL
      si.sort.cluster <- NULL
      for(q in 1:Q){
        si.sort <- c(si.sort,sort(si[cluster.list[[q]]],decreasing=T))
        si.sort.cluster <- c(si.sort.cluster,rep(q,length(cluster.list[[q]])))
      }
      return(list(cluster=si.sort.cluster,silhouette=si.sort,
                  silhouette.means=si.sort.cluster.means,silhouette.mean=si.mean))
    }
  }
  if(print.dims)if(is.null(A)){
    packageStartupMessage(
      sprintf("Y(%d,%d)~X(%d,%d)B(%d,%d)...",
              nrow(Y),ncol(Y),nrow(Y),Q,Q,ncol(Y)),appendLF=FALSE)
  }else{
    packageStartupMessage(
      sprintf("Y(%d,%d)~X(%d,%d)C(%d,%d)A(%d,%d)=XB(%d,%d)...",
              nrow(Y),ncol(Y),nrow(Y),Q,Q,nrow(A),nrow(A),ncol(Y),Q,ncol(Y)),appendLF=FALSE)
  }
  start.time <- Sys.time()
  set.seed(seed)
  if(is.vector(Y)) Y <- matrix(Y,nrow=1)
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  if(!is.null(A)){
    if(min(A)<0){
      warning("The matrix A should be non-negative.")
      stop()
    }
  }
  if(min(Y)<0){
    warning("The matrix Y should be non-negative.")
    stop()
  }
  is.X.scalar <- FALSE
  if(nrow(Y)>=2){
    if(min(nrow(Y),ncol(Y))>=Q){
      if(ncol(Y)==Q){
        X <- Y
      }else{
        res.kmeans <- stats::kmeans(t(Y),centers=Q,iter.max=maxit,nstart=nstart)
        X <- t(res.kmeans$centers)
      }
    }else{
      warning("It does not hold Q<=min(P,N) where dim(Y)=(P,N).")
      stop()
    }
  }else{
    X <- matrix(data=1,nrow=1,ncol=1)
    is.X.scalar <- TRUE
  }
  if(!is.X.scalar){
    if(X.restriction=="colSums") X <- t(t(X)/colSums(X))
    if(X.restriction=="colSqSums") X <- t(t(X)/colSums(X^2)^0.5)
    if(X.restriction=="totalSum") X <- X/sum(X)
  }
  if(is.null(A)) C <- matrix(1,nrow=ncol(X),ncol=ncol(Y)) else C <- matrix(1,nrow=ncol(X),ncol=nrow(A))
  objfunc.iter <- 0*(1:maxit)
  for(i in 1:maxit){
    if(is.null(A)) B <- C else B <- C %*% A
    XB <- X %*% B
    if(print.trace&i %% 10==0) print(paste0(format(Sys.time(), "%X")," ",i,"..."))
    if(method=="EU"){
      if(!is.X.scalar){
        X <- X*z((Y%*%t(B))/(XB%*%t(B)))
        if(X.restriction=="colSums") X <- t(t(X)/colSums(X))
        if(X.restriction=="colSqSums") X <- t(t(X)/colSums(X^2)^0.5)
        if(X.restriction=="totalSum") X <- X/sum(X)
      }
      if(is.null(A)) C <- C*z(crossprod(X,Y)/(crossprod(X,XB)+gamma*C)) else
        C <- C*z(crossprod(X,tcrossprod(Y,A))/(crossprod(X,tcrossprod(XB,A))+gamma*C))
      #if(is.null(A)) C <- C*z((t(X)%*%Y)/(t(X)%*%XB+gamma*C)) else
      #  C <- C*z((t(X)%*%Y%*%t(A))/(t(X)%*%XB%*%t(A)+gamma*C))
      objfunc.iter[i] <- sum((Y-XB)^2)+gamma*sum(C^2)
    }else{
      if(!is.X.scalar){
        X <- t(t(X*z(Y/XB)%*%t(B))/rowSums(B))
        if(X.restriction=="colSums") X <- t(t(X)/colSums(X))
        if(X.restriction=="colSqSums") X <- t(t(X)/colSums(X^2)^0.5)
        if(X.restriction=="totalSum") X <- X/sum(X)
      }
      if(is.null(A)) C <- C*z(crossprod(X,z(Y/XB))/(colSums(X)%o%rep(1,ncol(Y))+2*gamma*C)) else
        C <- C*z(crossprod(X,tcrossprod(z(Y/XB),A))/(colSums(X)%o%rowSums(A)+2*gamma*C))
      #if(is.null(A)) C <- C*(t(X)%*%z(Y/XB)/(colSums(X)%o%rep(1,ncol(Y))+2*gamma*C)) else
      #  C <- C*(t(X)%*%z(Y/XB)%*%t(A)/(colSums(X)%o%rowSums(A)+2*gamma*C))
      objfunc.iter[i] <- sum(-Y*z(log(XB))+XB)+gamma*sum(C^2)
    }
    if(i>=10){
      epsilon.iter <- abs(objfunc.iter[i]-objfunc.iter[i-1])/(abs(objfunc.iter[i])+0.1)
      if(epsilon.iter <= abs(epsilon)){
        objfunc.iter <- objfunc.iter[10:i]
        break
      }
    }
  }
  if(is.null(A)) B <- C else B <- C %*% A
  XB <- X %*% B
  if(method=="EU"){
    objfunc <- sum((Y-XB)^2)+gamma*sum(C^2)
  }else{
    objfunc <- sum(-Y*z(log(XB))+XB)+gamma*sum(C^2)
  }
  if(ncol(X)>1 & sum(rowSums(X)==1)==nrow(X)){
    index <- order(matrix(1:nrow(X)/nrow(X),nrow=1) %*% X)
    X <- X[,index]
    B <- B[index,]
    C <- C[index,]
  }
  rownames(C) <- paste0(prefix,1:nrow(C))
  rownames(X) <- rownames(Y)
  colnames(X) <- paste0(prefix,1:ncol(X))
  rownames(B) <- paste0(prefix,1:nrow(B))
  colnames(B) <- colnames(Y)
  r2 <- stats::cor(as.vector(XB),as.vector(Y))^2
  if(method=="EU"){
    ICp <- log(objfunc/prod(dim(Y)))+Q*sum(dim(Y))/prod(dim(Y))*log(prod(dim(Y))/sum(dim(Y)))
  }else{
    ICp <- NA
  }
  if(save.memory==FALSE){
    B.prob <- t(z(t(B)/colSums(B)))
    B.prob.sd.min <- min(apply(B.prob,1,stats::sd))
    B.cluster <- apply(B.prob,2,which.max)
    B.cluster[colSums(B.prob)==0] <- NA
    rownames(B.prob) <- paste0(prefix,1:nrow(B))
    colnames(B.prob) <- colnames(Y)
    colnames(XB) <- colnames(Y)
    X.prob <- z(X/rowSums(X))
    X.cluster <- apply(X.prob,1,which.max)
    X.cluster[rowSums(X.prob)==0] <- NA
    if(save.time){
      silhouette <- NA
      CPCC <- NA
    }else{
      silhouette <- mysilhouette(B.prob,B.cluster)
      if(Q>=2){
        M <- t(B.prob) %*% B.prob
        h.dist <- as.matrix(stats::cophenetic(stats::hclust(stats::as.dist(1-M))))
        up <- upper.tri(M)
        CPCC <- stats::cor(h.dist[up],(1-M)[up])
      }else{
        CPCC <- NA
      }
    }
  }else{
    B.prob <- NA
    B.prob.sd.min <- NA
    B.cluster <- NA
    XB <- NA
    X.prob <- NA
    X.cluster <- NA
    silhouette <- NA
    CPCC <- NA
  }
  if(epsilon.iter > abs(epsilon)) warning(paste0(
    "maximum iterations (",maxit,
    ") reached and the optimization hasn't converged yet."))
  end.time <- Sys.time()
  diff.time <- difftime(end.time,start.time,units="sec")
  diff.time.st <- ifelse(diff.time<=180,paste0(round(diff.time,1),"sec"),
                         paste0(round(diff.time/60,1),"min"))
  if(print.dims) packageStartupMessage(diff.time.st)
  result <- list(X=X,B=B,XB=XB,C=C,
                 B.prob=B.prob,B.cluster=B.cluster,
                 X.prob=X.prob,X.cluster=X.cluster,
                 objfunc=objfunc,objfunc.iter=objfunc.iter,r.squared=r2,
                 criterion=list(B.prob.sd.min=B.prob.sd.min,ICp=ICp,silhouette=silhouette,CPCC=CPCC))
  class(result) <- "nmfkc"
  return(result)
}


#' @title plot for return value of nmfkc function
#' @description \code{plot.nmfkc} plot for return value of nmfkc function
#' @param x return value of nmfkc function
#' @param ... arguments to be passed to plot function.
#' @export
plot.nmfkc <- function(x,...){
  plot(x$objfunc.iter,xlab="iter",ylab="objfunc",main=paste0("r.squared=",round(x$r.squared,3)),...)
}


#' @title predict for return value of nmfkc function
#' @description \code{predict.nmfkc} predict for return value of nmfkc function
#' @param x return value of nmfkc function
#' @param newA optionally, a new covariate matrix.
#' @param type The default is "response" given by the product of matrices X and B.
#'  If type is "prob", B.prob is used instead of B.
#'  If type is "class", class to maximize columns in B.prob.
#' @export
predict.nmfkc <- function(x,newA=NULL,type="response"){
  z <- function(x){
    x[is.nan(x)] <- 0
    x[is.infinite(x)] <- 0
    return(x)
  }
  if(is.null(newA)){
    if(type=="response"){
      result <- x$X %*% x$B
    }else{
      XB.prob <- x$X %*% x$B.prob
      if(type=="prob"){
        result <- XB.prob
      }else{
        result <- rownames(x$X)[apply(XB.prob,2,which.max)]
      }
    }
  }else{
    B <- x$C %*% newA
    if(type=="response"){
      result <- x$X %*% B
    }else{
      B.prob <- t(z(t(B)/colSums(B)))
      XB.prob <- x$X %*% B.prob
      if(type=="prob"){
        result <- XB.prob
      }else{
        result <- rownames(x$X)[apply(XB.prob,2,which.max)]
      }
    }
  }
  return(result)
}


#' @title class matrix for categorical vector
#' @description \code{predict.nmfkc} predict for return value of nmfkc function
#' @param x categorical vector or vector of factor
#' @export
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' Y <- nmfkc.class(iris$Species)
#' Y[,1:6]
nmfkc.class <- function(x){
  if(!is.factor(x)) x <- as.factor(x)
  unix <- levels(x)
  X <- matrix(0,nrow=length(unix),ncol=length(x))
  rownames(X) <- unix
  for(j in 1:length(unix)) X[j,] <- ifelse(x==unix[j],1,0)
  result <- X
  return(result)
}

#' @title Normalizing to a value between 0 and 1 for matrix
#' @description \code{nmfkc.normalize} Normalize x using the minimum and maximum values of each column of the reference matrix
#' @param x matrix to be normalized
#' @param ref matrix referencing minimum and maximum values.
#' @export
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' x <- nmfkc.normalize(iris[,-5])
#' apply(x,2,range)
nmfkc.normalize <- function(x,ref=x){
  if(is.vector(x)==TRUE){
    x <- matrix(x,ncol=1)
    ref <- matrix(ref,ncol=1)
  }
  r <- apply(ref,2,range)
  y <- 0*x
  for(j in 1:ncol(x))y[,j] <- (x[,j]-r[1,j])/(r[2,j]-r[1,j])
  return(y)
}


#' @title Performing k-fold cross validation on NMF (Non-negative Matrix Factorization) with Kernel Covariate
#' @description \code{nmfkc.cv} apply cross validation method for k-partitioned columns of Y~XCA=XB
#'  where observation matrix Y(P,N),
#'  covariate matrix A(R,N),
#'  basis matrix X(P,Q) and Q<=min(P,N),
#'  parameter matrix C(Q,R)
#'  and coefficient matrix B(Q,N).
#'  Note that Y(N,P) and A(R,N) are given, and X(P,Q) and C(Q,R) are unknown.
#' @param Y observation matrix
#' @param A covariate matrix. Without covariate, identity matrix is used.
#' @param Q rank of basis matrix and Q<=min(P,N) where Y(P,N)
#' @param div number of partition usually described as "k" of k-fold
#' @param seed integer used as argument in set.seed function
#' @param ... arguments to be passed to nmfkc function.
#'  which controls the reproducibility of the partition.
#' @return objfunc: last objective function
#' @return objfunc.block: objective function at each block
#' @return block: partition block index (1,...,div) assigned to each column of Y
#' @export
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example 1.
#' Y <- matrix(cars$dist,nrow=1)
#' A <- rbind(1,cars$speed)
#' result <- nmfkc.cv(Y,A,Q=1)
#' result$objfunc
#'
#' # Example 2.
#' Y <- matrix(cars$dist,nrow=1)
#' U <- matrix(c(5,10,15,20,25),nrow=1)
#' V <- matrix(cars$speed,nrow=1)
#' betas <- 25:35/1000
#' objfuncs <- 0*(1:length(betas))
#' for(i in 1:length(betas)){
#'   print(i)
#'   A <- nmfkc.kernel(U,V,beta=betas[i])
#'   result <- nmfkc.cv(Y,A,Q=1,div=10)
#'   objfuncs[i] <- result$objfunc
#' }
#' min(objfuncs)
#' (beta.best <- betas[which.min(objfuncs)])
#' # objective function by beta
#' plot(betas,objfuncs,type="o",log="x")
#' table(result$block) # partition block of cv
#' # fitted curve
#' A <- nmfkc.kernel(U,V,beta=beta.best)
#' result <- nmfkc(Y,A,Q=1)
#' plot(as.vector(V),as.vector(Y))
#' lines(as.vector(V),as.vector(result$XB),col=2,lwd=2)

nmfkc.cv <- function(Y,A=NULL,Q=2,div=5,seed=123,...){
  arglist=list(...)
  gamma <- ifelse("gamma" %in% names(arglist),arglist$gamma,0)
  epsilon <- ifelse("epsilon" %in% names(arglist),arglist$epsilon,1e-4)
  maxit <- ifelse("maxit" %in% names(arglist),arglist$maxit,5000)
  method <- ifelse("method" %in% names(arglist),arglist$method,"EU")
  X.restriction <- ifelse("X.restriction" %in% names(arglist),arglist$X.restriction,"colSums")
  nstart <- ifelse("nstart" %in% names(arglist),arglist$nstart,1)
  prefix <- ifelse("prefix" %in% names(arglist),arglist$prefix,"Basis")
  print.trace <- ifelse("print.trace" %in% names(arglist),arglist$print.trace,FALSE)
  print.dims <- ifelse("print.dims" %in% names(arglist),arglist$print.dims,FALSE)
  save.time <- TRUE
  save.memory <- TRUE
  if(is.null(A)) A <- diag(ncol(Y))
  is.symmetric.matrix <- function(A){
    result <- FALSE
    if(nrow(A)==ncol(A)){
      result <- sum(abs(t(A)-A))==0
    }
    return(result)
  }
  is.identity.matrix <- function(A){
    result <- FALSE
    if(nrow(A)==ncol(A)&min(A)==0&max(A)==1){
      if(prod(diag(A))==1&sum(A-diag(nrow(A)))==0) result <- TRUE
    }
    return(result)
  }
  z <- function(x){
    x[is.nan(x)] <- 0
    x[is.infinite(x)] <- 0
    return(x)
  }
  optimize.B.from.Y <- function(result,Y,gamma,epsilon,maxit,method){
    X <- result$X
    C <- matrix(1,nrow=ncol(X),ncol=ncol(Y))
    oldSum <- 0
    for(l in 1:maxit){
      B <- C
      XB <- X %*% B
      if(method=="EU"){
        C <- C*z((t(X)%*%Y)/(t(X)%*%XB+gamma*C))
      }else{
        C <- C*(t(X)%*%z(Y/XB)/(colSums(X)%o%rep(1,ncol(Y))+2*gamma*C))
      }
      newSum <- sum(C)
      if(l>=10){
        epsilon.iter <- abs(newSum-oldSum)/(abs(newSum)+0.1)
        if(epsilon.iter <= abs(epsilon)) break
      }
      oldSum <- sum(C)
    }
    B <- C
    XB <- X %*% B
    colnames(B) <- colnames(Y)
    colnames(XB) <- colnames(Y)
    if(epsilon.iter > abs(epsilon)) warning(paste0(
      "maximum iterations (",maxit,
      ") reached and the optimization hasn't converged yet."))
    return(list(B=B,XB=XB))
  }
  is.identity <- is.identity.matrix(A)
  is.symmetric.matrix <- is.symmetric.matrix(A)
  n <- ncol(Y)
  (remainder <- n %% div)
  (division <- n %/% div)
  block <- 0*(1:n)
  set.seed(seed)
  index <- sample(1:n,n,replace=F)
  for(i in 1:(div-1)){
    plus <- ifelse(i<=remainder,1,0)
    j <- index[1:(division+plus)]
    block[j] <- i
    index <- index[-(1:(division+plus))]
  }
  block[index] <- div
  index <- block
  objfunc.block <- 0*(1:div)
  for(j in 1:div){
    Y_j <- Y[,index!=j]
    Yj <- Y[,index==j]
    if(is.symmetric.matrix){
      A_j <- A[index!=j,index!=j] # kernel matrix or identity matrix
    }else{
      A_j <- A[,index!=j] # ordinary design matrix
    }
    res_j <- nmfkc(Y_j,A_j,Q,gamma,epsilon,maxit,method,X.restriction,nstart,seed,prefix,print.trace,print.dims,save.time,save.memory)
    if(is.identity){
      resj <- optimize.B.from.Y(res_j,Yj,gamma,epsilon,maxit,method)
      XBj <- resj$XB
    }else{
      if(is.symmetric.matrix){
        Aj <- A[index!=j,index==j]
      }else{
        Aj <- A[,index==j]
      }
      XBj <- res_j$X %*% res_j$C %*% Aj
    }
    if(method=="EU"){
      objfunc.block[j] <- sum((Yj-XBj)^2)+gamma*sum(res_j$C^2)
    }else{
      objfunc.block[j] <- sum(-Yj*z(log(XBj))+XBj)+gamma*sum(res_j$C^2)
    }
  }
  objfunc <- sum(objfunc.block)
  return(list(objfunc=objfunc,objfunc.block=objfunc.block,block=index))
}


#' @title Rank selection diagnostics using figure
#' @description \code{nmfkc.rank} provides rank selection diagnostics. The method is still under development.
#' @param Y observation matrix
#' @param A covariate matrix. Without covariate, identity matrix is used.
#' @param Q vector of ranks to be diagnosed.
#' @param plot The default is plot=TRUE and draw a graph.
#' @param ... arguments to be passed to nmfkc function.
#' @export
#' @references Brunet, J.P., Tamayo, P., Golub, T.R., Mesirov, J.P. (2004) Metagenes and molecular pattern discovery using matrix factorization. Proc. Natl. Acad. Sci. USA 2004, 101, 4164–4169. \url{https://doi.org/10.1073/pnas.0308531101}
#' @references Punera, K. and Ghosh, J. (2008). CONSENSUS-BASED ENSEMBLES OF SOFT CLUSTERINGS. Applied Artificial Intelligence, 22(7–8), 780–810. \url{https://doi.org/10.1080/08839510802170546}
#' @examples
#' # install.packages("remotes")
#' # remotes::install_github("ksatohds/nmfkc")
#' # Example.
#' library(nmfkc)
#' Y <- t(iris[,-5])
#' nmfkc.rank(Y,Q=2:4)

nmfkc.rank <- function(Y,A=NULL,Q=2:min(5,ncol(Y),nrow(Y)),plot=TRUE,...){
  arglist=list(...)
  AdjustedRandIndex <- function(x){
    choose2 <- function(n) choose(n,2)
    a <- sum(apply(x,c(1,2),choose2))
    ab <- sum(sapply(rowSums(x),choose2))
    b <- ab-a
    ac <- sum(sapply(colSums(x),choose2))
    c <- ac-a
    total <- choose2(sum(x))
    d <- total-a-b-c
    (ri <- (a+d)/total)
    e <- ab*ac/total+(total-ab)*(total-ac)/total
    (ari <- (a+d-e)/(total-e))
    return(list(RI=ri,ARI=ari))
  }
  gamma <- ifelse("gamma" %in% names(arglist),arglist$gamma,0)
  epsilon <- ifelse("epsilon" %in% names(arglist),arglist$epsilon,1e-4)
  maxit <- ifelse("maxit" %in% names(arglist),arglist$maxit,5000)
  method <- ifelse("method" %in% names(arglist),arglist$method,"EU")
  X.restriction <- ifelse("X.restriction" %in% names(arglist),arglist$X.restriction,"colSums")
  seed <- ifelse("seed" %in% names(arglist),arglist$nstart,123)
  nstart <- ifelse("nstart" %in% names(arglist),arglist$nstart,1)
  prefix <- ifelse("prefix" %in% names(arglist),arglist$prefix,"Basis")
  print.trace <- ifelse("print.trace" %in% names(arglist),arglist$print.trace,FALSE)
  print.dims <- ifelse("print.dims" %in% names(arglist),arglist$print.dims,TRUE)
  save.time <- ifelse("save.time" %in% names(arglist),arglist$save.time,TRUE)
  save.memory <- FALSE
  r.squared <- 0*Q; names(r.squared) <- Q
  ICp <- 0*Q; names(ICp) <- Q
  silhouette <- 0*Q; names(silhouette) <- Q
  CPCC <- 0*Q; names(CPCC) <- Q
  B.prob.sd.min <- 0*Q; names(B.prob.sd.min) <- Q
  ARI <- 0*Q; names(ARI) <- Q
  for(q in 1:length(Q)){
    if(save.time){
      result <- nmfkc(Y,A,Q=Q[q],gamma,epsilon,maxit,method,X.restriction,nstart,seed,prefix,print.trace,print.dims,save.time=T,save.memory)
      CPCC[q] <- NA
      silhouette[q] <- NA
    }else{
      result <- nmfkc(Y,A,Q=Q[q],gamma,epsilon,maxit,method,X.restriction,nstart,seed,prefix,print.trace,print.dims,save.time=F,save.memory)
      CPCC[q] <- result$criterion$CPCC
      silhouette[q] <- result$criterion$silhouette$silhouette.mean
    }
    r.squared[q] <- result$r.squared
    ICp[q] <- result$criterion$ICp
    if(q==1){
      ARI[q] <- NA
      cluster.old <- result$B.cluster
    }else{
      df <- data.frame(old=cluster.old,new=result$B.cluster)
      df <- df[stats::complete.cases(df),]
      f <- table(df$old,df$new)
      ARI[q] <- AdjustedRandIndex(f)$ARI
      cluster.old <- result$B.cluster
    }
    B.prob.sd.min[q] <- result$criterion$B.prob.sd.min
  }
  if(plot){
    # Criterion
    plot(Q,r.squared,type="l",col=2,xlab="Rank",ylab="Criterion",ylim=c(0,1),lwd=3)
    graphics::text(Q,r.squared,Q)
    legend <- "r.squared"
    fill <- 2
    graphics::lines(Q,B.prob.sd.min,col=3,lwd=3)
    graphics::text(Q,B.prob.sd.min,Q)
    legend <- c(legend,"B.prob.sd.min")
    fill <- c(fill,3)
    graphics::lines(Q[-1],ARI[-1],col=4,lwd=3)
    graphics::text(Q[-1],ARI[-1],Q[-1])
    legend <- c(legend,"ARI for Q-1")
    fill <- c(fill,4)
    if(!save.time){
      graphics::lines(Q,silhouette,col=7,lwd=3)
      graphics::text(Q,silhouette,Q)
      legend <- c(legend,"silhouette")
      fill <- c(fill,7)
      graphics::lines(Q,CPCC,col=6,lwd=3)
      graphics::text(Q,CPCC,Q)
      legend <- c(legend,"CPCC")
      fill <- c(fill,6)
    }
    graphics::legend("bottomleft",legend=legend,fill=fill,bg=NULL)
  }
  invisible(data.frame(Q,r.squared,ICp,B.prob.sd.min,ARI,silhouette,CPCC))
}

