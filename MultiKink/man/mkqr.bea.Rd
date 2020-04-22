\name{mkqr.bea}
\alias{mkqr.bea}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit the multi-kink quantile regression in absence of the number of change points.}
\description{A quantile BIC combining the backward elimination algorithm to estimate the multi-kink quantile regression.}
\usage{
mkqr.bea(y,thre.x,cont.z,tau,Cn,bandwidth_type="Hall-Sheather",control=fit.control())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{A vector of response.}
  \item{thre.x}{A scalar covariate with threshold effect.}
  \item{cont.z}{A vector of covariates with constant slopes.}
  \item{tau}{The quantile level that belongs to (0,1).}
  \item{Cn}{A psotive number corresponding to different types of BIC.}
  \item{bandwidth_type}{The bandwidth type. Specify one from "Hall-Sheather","Bofinger","Chamberlain". Default is "Hall-Sheather".}
  \item{control}{A list returned by \code{fit.control}.}
}
\details{
this function does not need to set the number of change points in advance. The number of change points will be estimated automatically.
}
\value{
A list with the elements
\item{bet.est}{The estimated regression coefficients with intercept.}
\item{bet.se}{The estimated standard error of the regression coefficients.}
\item{psi.est}{The estimated change points.}
\item{psi.se}{The estiamted standard errors of threshold parameters.}
\item{n.psi}{The estiamted number of change points.}
}
\references{
Wei Zhong, Chuang Wan and Wenyang Zhang. (2020) Estimation and inference for multi-kink quantile regression. working paper.
}
\seealso{
\code{\link{mkqr.fit}}, \code{\link{fit.control}}
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
fit <- mkqr.bea(y,x,z,tau,Cn=log(n))


## The example of triceps dataset
data(triceps)
y <- triceps$lntriceps
x <- triceps$age
n <- length(y)
z <- NULL
tau <- 0.5
fit <- mkqr.bea(y, x, z, tau,Cn=log(n),control=fit.control(K.max=4))
proc.time() - ptm
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mkqr.bea }
