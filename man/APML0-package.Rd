\name{APML0-package}
\alias{APML0-package}
\docType{package}
\title{
Augmented and Penalized Minimization Method L0
}
\description{
Fit linear and Cox models regularized with L0, lasso (L1), elastic-net (L1 and L2), or net (L1 and Laplacian) penalty, and their adaptive forms, such as adaptive lasso / elastic-net and net adjusting for signs of linked coefficients. It solves L0 penalty problem by simultaneously selecting regularization parameters and the number of non-zero coefficients. This augmented and penalized minimization method provides an approximation solution to the L0 penalty problem, but runs as fast as L1 regularization problem.

The package uses one-step coordinate descent algorithm and runs extremely fast by taking into account the sparsity structure of coefficients. It could deal with very high dimensional data.
}
\details{
\tabular{ll}{
Package: \tab APML0\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2017-07-15\cr
License: \tab GPL (>= 2)\cr
}
Functions:
\code{\link{APML0}}, \code{\link{print.APML0}}\cr
}
\author{
Xiang Li, Shanghong Xie, Donglin Zeng and Yuanjia Wang\cr
Maintainer: Xiang Li <xli256@its.jnj.com>
}
\references{Boyd, S., Parikh, N., Chu, E., Peleato, B., & Eckstein, J. (2011).
  \emph{Distributed optimization and statistical learning via the alternating direction method of multipliers. Foundations and Trends in Machine Learning, 3(1), 1-122.}\cr
  \url{http://dl.acm.org/citation.cfm?id=2185816}\cr
  Friedman, J., Hastie, T. and Tibshirani, R. (2010).
  \emph{Regularization paths for generalized linear models via coordinate descent, Journal of Statistical Software, Vol. 33(1), 1.}\cr
  \url{http://www.jstatsoft.org/v33/i01/}\cr
  Li, C., and Li, H. (2010).
  \emph{Variable selection and regression analysis for graph-structured covariates with an application to genomics. The annals of applied statistics, 4(3), 1498.}\cr
  \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3423227/}\cr
   Sun, H., Lin, W., Feng, R., and Li, H. (2014)
  \emph{Network-regularized high-dimensional cox regression for analysis of genomic data, Statistica Sinica.}\cr
  \url{http://www3.stat.sinica.edu.tw/statistica/j24n3/j24n319/j24n319.html}\cr
}

\keyword{Number of non-zeros}
\keyword{Regularization}
\keyword{Package}

\examples{
###  Linear model  ###
set.seed(1213)
N=100;p=30;p1=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(p1)
xb=x[,1:p1]\%*\%beta
y=rnorm(N,xb)

fiti=APML0(x,y,penalty="Lasso",nlambda=10) # Lasso
fiti2=APML0(x,y,penalty="Lasso",nlambda=10,nfolds=10) # Lasso
# attributes(fiti)


###  Cox model  ###
set.seed(1213)
N=100;p=30;p1=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(p1)
xb=x[,1:p1]\%*\%beta
ty=rexp(N,exp(xb))
tcens=rbinom(n=N,prob=.3,size=1)  # censoring indicator
y=cbind(time=ty,status=1-tcens)

fiti=APML0(x,y,family="cox",penalty="Lasso",nlambda=10) # Lasso
fiti2=APML0(x,y,family="cox",penalty="Lasso",nlambda=10,nfolds=10) # Lasso
# attributes(fiti)
}

