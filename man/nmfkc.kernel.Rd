% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nmfkc.R
\name{nmfkc.kernel}
\alias{nmfkc.kernel}
\title{Creating kernel matrix from covariates}
\source{
Satoh, K. (2024) Applying Non-negative Matrix Factorization with Covariates to the Longitudinal Data as Growth Curve Model. arXiv preprint arXiv:2403.05359. \url{https://arxiv.org/abs/2403.05359}
}
\usage{
nmfkc.kernel(U, V = U, beta)
}
\arguments{
\item{U}{covariate matrix U(K,N)=(u_1,...,u_N) each row might be normalized in advance}

\item{V}{covariate matrix V(K,M)=(v_1,...,v_M) usually used for prediction, and the default value is U.}

\item{beta}{parameter of kernel matrix A of which element is defined by Gauss kernel function, exp(-beta*|u_n-v_m|^2)}
}
\value{
kernel matrix A(N,M)
}
\description{
\code{nmfkc.kernel} create kernel matrix from covariate matrix
}
\examples{
# install.packages("remotes")
# remotes::install_github("ksatohds/nmfkc")
# Example 1.
U <- matrix(1:3,nrow=1,ncol=3)
print(U)
A <- nmfkc.kernel(U,beta=1)
print(A)
print(log(A))

# Example 2.
Y <- matrix(cars$dist,nrow=1)
U <- matrix(c(5,10,15,20,25),nrow=1)
V <- matrix(cars$speed,nrow=1)
A <- nmfkc.kernel(U,V,beta=0.031) # see, nmfkc.cv
result <- nmfkc(Y,A,Q=1)
plot(as.vector(V),as.vector(Y))
lines(as.vector(V),as.vector(result$XB),col=2,lwd=2)
}
\references{
Ding, C., Li, T., Peng, W. and Park, H. (2006) Orthogonal Nonnegative Matrix Tri-Factorizations for Clustering, Proceedings of the 12th ACM SIGKDD international conference on Knowledge discovery and data mining, 126-135. \url{https://doi.org/10.1145/1150402.1150420}

Potthoff, R.F., and Roy, S.N. (1964). A generalized multivariate analysis of variance model useful especially for growth curve problems. Biometrika, 51, 313-326. \url{https://doi.org/10.2307/2334137}
}
