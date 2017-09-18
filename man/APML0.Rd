\name{APML0}
\alias{APML0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a Model with Various Regularization Forms
}
\description{
Fit linear, logistic and Cox models regularized with L0, lasso (L1), elastic-net (L1 and L2), or net (L1 and Laplacian) penalty, and their adaptive forms, such as adaptive lasso / elastic-net and net adjusting for signs of linked coefficients. It solves L0 penalty problem by simultaneously selecting regularization parameters and the number of non-zero coefficients. This augmented and penalized minimization method provides an approximation solution to the L0 penalty problem, but runs as fast as L1 regularization problem.

The package uses one-step coordinate descent algorithm and runs extremely fast by taking into account the sparsity structure of coefficients. It could deal with very high dimensional data.
}
\usage{
APML0(x, y, family=c("gaussian", "binary", "cox"), penalty=c("Lasso","Enet", "Net"),
Omega=NULL, alpha=1.0, lambda=NULL, nlambda=50, rlambda=NULL,
wbeta=rep(1,ncol(x)), sgn=rep(1,ncol(x)), nfolds=1, foldid=NULL, inzero=TRUE,
isd=FALSE, iysd=FALSE, keep.beta=FALSE, ifast=TRUE, thresh=1e-7, maxit=1e+5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{input matrix. Each row is an observation vector.
}
  \item{y}{response variable. For \code{family = "gaussian"},
  \code{y} is a continuous vector. For \code{family = "binary"},
  \code{y} is a binary vector with 0 and 1. For \code{family = "cox"}, \code{y} is a two-column matrix with columns named `time' and `status'. `status' is a binary variable, with `1' indicating event, and `0' indicating right censored.
}
  \item{family}{type of outcome. Can be "gaussian", "binary" or "cox".
}
  \item{penalty}{penalty type. Can choose \code{"Net"}, \code{"Enet"} (elastic net) and \code{"Lasso"}. For \code{"Net"}, need to specify \code{Omega}; otherwise, \code{"Enet"} is performed. For \code{penalty = "Net"}, the penalty is defined as \deqn{\lambda*{\alpha*||\beta||_1+(1-\alpha)/2*(\beta^{T}L\beta)},}
where \eqn{L} is a Laplacian matrix calculated from \code{Omega}.
}
  \item{Omega}{adjacency matrix with zero diagonal and non-negative off-diagonal, used for \code{penalty = "Net"} to calculate Laplacian matrix.
}
  \item{alpha}{ratio between L1 and Laplacian for \code{"Net"}, or between L1 and L2 for \code{"Enet"}. Default is \code{alpha = 1.0}, i.e. lasso.
}
  \item{lambda}{a user supplied decreasing sequence. If \code{lambda = NULL}, a sequence of \code{lambda} is generated based on \code{nlambda} and \code{rlambda}. Supplying a value of \code{lambda} overrides this.
}
  \item{nlambda}{number of \code{lambda} values. Default is 50.
}
  \item{rlambda}{fraction of \code{lambda.max} to determine the smallest value for \code{lambda}. The default is \code{rlambda = 0.0001} when the number of observations is larger than or equal to the number of variables; otherwise, \code{rlambda = 0.01}.
}
  \item{wbeta}{penalty weights used with L1 penalty (adaptive L1), given by \eqn{\sum_{j=1}^qw_j|\beta_j|.} The \code{wbeta} is a vector of non-negative values and works as adaptive L1. No penalty is imposed for those coefficients with zero values in \code{wbeta}. Default is 1 for all coefficients.
}
  \item{sgn}{sign adjustment used with Laplacian penalty (adaptive Laplacian). The \code{sgn} is a vector of 1 or -1. The \code{sgn} could be based on an initial estimate of \eqn{\beta}, and 1 is used for \eqn{\beta>0} and -1 is for \eqn{\beta<0}. Default is 1 for all coefficients.
}
  \item{nfolds}{number of folds. With \code{nfolds = 1} and \code{foldid = NULL} by default, cross-validation is not performed. For cross-validation, smallest value allowable is \code{nfolds = 3}. Specifying \code{foldid} overrides \code{nfolds}.
}
  \item{foldid}{an optional vector of values between 1 and \code{nfolds} specifying which fold each observation is in.
}
  \item{inzero}{logical flag for simultaneously selecting the number of non-zero coefficients with \code{lambda}. Default is \code{inzero = TRUE}.
}
  \item{isd}{logical flag for outputting standardized coefficients. \code{x} is always standardized prior to fitting the model. Default is \code{isd = FALSE}, returning \eqn{\beta} on the original scale.
}
  \item{iysd}{logical flag for standardizing \code{y} prior to computation, for \code{family = "gaussian"}. The returning coefficients are always based the original \code{y} (unstandardized). Default is \code{isd = FALSE}.
}
  \item{keep.beta}{logical flag for returning estimates for all \code{lambda} values. For \code{keep.beta = FALSE}, only return the estimate with the minimum cross-validation value.
}
  \item{ifast}{logical flag for efficient calculation of risk set updates for \code{family = "cox"}. Default is \code{ifast = TRUE}.
}
  \item{thresh}{convergence threshold for coordinate descent. Default value is \code{1E-7}.
}
  \item{maxit}{Maximum number of iterations for coordinate descent. Default is \code{10^5}.
}
}
\details{
  One-step coordinate descent algorithm is applied for each \code{lambda}. For \code{family = "cox"}, \code{ifast = TRUE} adopts an efficient way to update risk set and sometimes the algorithm ends before all \code{nlambda} values of \code{lambda} have been evaluated. To evaluate small values of \code{lambda}, use \code{ifast = FALSE}. The two methods only affect the efficiency of algorithm, not the estimates.

  \code{x} is always standardized prior to fitting the model and the estimate is returned on the original scale. For \code{family = "gaussian"}, y is centered by removing its mean, so there is no intercept output.

  Cross-validation is used for tuning parameters. For \code{inzero =  TRUE}, we further select the number of non-zero coefficients obtained from regularized model at each \code{lambda}. This is motivated by formulating L0 variable selection in an augmented form, which shows significant improvement over the commonly used regularized methods without this technique.
}
\value{
An object with S3 class \code{"APML0"}.
  \item{Beta}{a sparse Matrix of coefficients, stored in class "dgCMatrix".}
  \item{Beta0}{coefficients after additionally tuning the number of non-zeros, for \code{inzero = TRUE}.}
  \item{fit}{a data.frame containing \code{lambda} and the number of non-zero coefficients \code{nzero}. With cross-validation, additional results are reported, such as average cross-validation partial likelihood \code{cvm} and its standard error \code{cvse}, and \code{index} with `*' indicating the minimum \code{cvm}. For \code{family = "gaussian"}, \code{rsq} is also reported.}
  \item{fit0}{a data.frame containing \code{lambda}, \code{cvm} and \code{nzero} based on \code{inzero = TRUE}. \code{cvm} in \code{fit0} may be different from \code{cvm} in \code{fit}, because the constaint on the number of non-zeros is imposed in the cross-validation. The maximum number of non-zeros is based on the full dataset not the one used in the cross-validation.}

  \item{lambda.min}{value of \code{lambda} that gives minimum \code{cvm}.}
  \item{lambda.opt}{value of \code{lambda} based on \code{inzero = TRUE}.}

  \item{penalty}{penalty type.}
  \item{adaptive}{logical flags for adaptive version (see above).}
  \item{flag}{convergence flag (for internal debugging). \code{flag = 0} means converged.}
}
\references{Li, X., Xie, S., Zeng, D., Wang, Y. (2017).
\emph{Efficient Method to Optimally Identify Important Biomarkers for Disease Outcomes with High-dimensional Data. Statistics in Medicine, accepted.}\cr
Boyd, S., Parikh, N., Chu, E., Peleato, B., & Eckstein, J. (2011).
  \emph{Distributed optimization and statistical learning via the alternating direction method of multipliers. Foundations and Trends in Machine Learning, 3(1), 1-122.}\cr
  \url{http://dl.acm.org/citation.cfm?id=2185816}\cr
  Friedman, J., Hastie, T. and Tibshirani, R. (2010).
  \emph{Regularization paths for generalized linear models via coordinate descent, Journal of Statistical Software, Vol. 33(1), 1.}\cr
  \url{http://www.jstatsoft.org/v33/i01/}\cr
}

\author{
Xiang Li, Shanghong Xie, Donglin Zeng and Yuanjia Wang\cr
Maintainer: Xiang Li <xli256@its.jnj.com>, Shanghong Xie <sx2168@cumc.columbia.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\section{Warning}{
  It may terminate and return \code{NULL}.
}

\seealso{
\code{\link{print.APML0}}
}
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



###  Logistic model  ###
set.seed(1213)
N=100;p=30;p1=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(p1)
xb=x[,1:p1]\%*\%beta
y=rbinom(n=N, size=1, prob=1.0/(1.0+exp(-xb)))

fiti=APML0(x,y,family="binary",penalty="Lasso",nlambda=10) # Lasso
fiti2=APML0(x,y,family="binary",penalty="Lasso",nlambda=10,nfolds=10) # Lasso
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Number of non-zeros}
\keyword{Regularization}


