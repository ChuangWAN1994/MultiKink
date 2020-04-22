\name{kinkTest}
\alias{kinkTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test the existence of kink effect in the multi-kink quantile regression}
\description{This function is used to  calculate the score test statistics and its asymptotical p-value is obtained by using wild bootstrap.}
\usage{
kinkTest(y,thre.x,cont.z,tau,NB=200,sparsity="nid",
                     bandwidth_type="Hall-Sheather")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{A vector of response.}
  \item{thre.x}{A scalar covariate with threshold effect.}
  \item{cont.z}{A vector of covariates with constant slopes.}
  \item{tau}{A given quantile level belongs to (0,1).}
  \item{NB}{Resampling times, 200 for default.}
  \item{sparsity}{The error term type. Specify one from "iid" and "nid". Default is "nid".}
  \item{bandwidth_type}{The bandwidth type. Specify one from "Hall-Sheather","Bofinger","Chamberlain". Default is "Hall-Sheather".}
}

\value{
A list with the elements
\item{Tn}{The statistic based on original data.}
\item{Tn.NB}{The statistics by wild bootstrap.}
\item{pv}{The p-value by wild bootstrap.}
}
\references{
Wei Zhong, Chuang Wan and Wenyang Zhang. (2020) Estimation and inference for multi-kink quantile regression. working paper.
}
\author{
Chuang Wan
}

\examples{
ptm <- proc.time()
##Simulated data
set.seed(123)
n <- 500
tau <- 0.5
x <- runif(n,-5,5)
z <- rnorm(n,1,1)
e <- rnorm(n,0,1)-qnorm(tau,0,1)
y <- 1+x-0.5*pmax(x-2,0)+z+e
obj <- kinkTest(y, x, z, tau, NB = 200,sparsity="nid",bandwidth_type="Hall-Sheather")
proc.time() - ptm
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ kinkTest }% use one of  RShowDoc("KEYWORDS")
