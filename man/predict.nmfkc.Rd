% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nmfkc.R
\name{predict.nmfkc}
\alias{predict.nmfkc}
\title{predict for return value of nmfkc function}
\usage{
\method{predict}{nmfkc}(x, newA = NULL, type = "response")
}
\arguments{
\item{x}{return value of nmfkc function}

\item{newA}{optionally, a new covariate matrix.}

\item{type}{The default is "response" given by the product of matrices X and B.
If type is "prob", B.prob is used instead of B.
If type is "class", class to maximize columns in B.prob.}
}
\description{
\code{predict.nmfkc} predict for return value of nmfkc function
}
