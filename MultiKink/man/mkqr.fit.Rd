\name{mkqr.fit}
\alias{mkqr.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit the multi-kink quantile regression conditonal on a given or pre-specified number of change points.
}
\description{
The bootstrap restarting iterative segmented quantile algorithm for the multi-kink quantile regression.
}
\usage{
mkqr.fit(y,thre.x,cont.z,psi,k,tau,bandwidth_type="Hall-Sheather",control=fit.control())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{A vector of response.}
  \item{thre.x}{A scalar covariate with threshold effect.}
  \item{cont.z}{A vector of covariates with constant slopes.}
  \item{psi}{numeric vector to indicate the starting values for the changepoints. When \code{psi=NULL} (default), \code{k} quantiles are assumed.}
  \item{k}{The pre-specified number of change points.}
  \item{tau}{The quantile level that belongs to (0,1).}
  \item{bandwidth_type}{The bandwidth type. Specify one from "Hall-Sheather","Bofinger","Chamberlain". Default is "Hall-Sheather".}
  \item{control}{A list returned by \code{fit.control}.}
}
\details{either 'psi' or 'k' must be given.}
\value{
A list with the elements
\item{bet.est}{The estimated regression coefficients with intercept.}
\item{bet.se}{The estimated standard error of the regression coefficients.}
\item{psi.est}{The estimated change points.}
\item{psi.se}{The estiamted standard errors of threshold parameters.}
}
\references{
Wood S.N. (2001) Minimizing model fitting objectives that contain spurious local minima by bootstrap restarting. Biometrics, 57(1): 240-244.

Muggeo, V.M.R., Adelfio, G. (2011) Efficient change point detection in genomic sequences of continuous measurements. Bioinformatics 27, 161-166.

Wei Zhong, Chuang Wan and Wenyang Zhang. (2020) Estimation and inference for multi-kink quantile regression. working paper.
}
\seealso{
\code{\link{mkqr.bea}},\code{\link{fit.control}}
}
\author{
Chuang Wan
}


\examples{
\dontrun{
ptm <- proc.time()
##simulated data
set.seed(123)
n <- 500
tau <- 0.5
bet0 <- c(1,1,1,-3,4,4)
psi0 <- c(-3,0,3)
k <- length(psi0)
x <- runif(n, -5, 5)
z <- rnorm(n, 1, 1)
X <- matrix(rep(x,k),nrow=n)
PSI <- matrix(rep(psi0,rep(n,k)),ncol=k)
XZ <- cbind(1,z,x,pmax((X-PSI),0))
e <- rnorm(n,0,1) - qnorm(tau, 0,1)
y <- as.vector(XZ \%*\% bet0) + e
fit <- mkqr.fit(y,x,z,psi=NULL,k=3,tau)

## The example of triceps dataset
data(triceps)
y <- triceps$lntriceps
x <- triceps$age
n <- length(y)
z <- NULL
tau <- 0.5
fit <- mkqr.fit(y, x, z, psi=NULL,k=2,tau)
proc.time() - ptm
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mkqr.fit }
