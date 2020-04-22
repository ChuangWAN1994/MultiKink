\name{fit.control}
\alias{fit.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auxiliary parameters to  control the model fitting.
}
\description{
This is a  user interface to control the auxiliary parameter in fitting multi-kink quantile regression.
}
\usage{
fit.control(toll=1e-4,h=1,it.max=50,K.max=6,stop.if.error=TRUE,dev0=NULL,visual=FALSE,
           visualBoot=FALSE,pow=c(1,1),digits=NULL,grid=NULL,n.boot=20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{toll}{Positive convergence tolerance.}
  \item{h}{Positive factor (from zero to one) modifying the increments in kink parameter updates during  the iterative process.}
  \item{it.max}{Positive integer for the maximal number of iterations.}
  \item{K.max}{Positive integer for the maximal given number of kink points.}
  \item{stop.if.error}{Logical indicating if the estimation algorithm shoud be stopped if some kink point estimators belong to the non-admissible set. Default is \code{FALSE} which suggests to remove the non-addimissble change points automatically.}
  \item{dev0}{Initial objective value or deviance. Default is \code{NULL} which implies that the initial values is unknown.}
  \item{visual}{Logical indicating if the results of  estimation process should be printed at each iteration.}
  \item{visualBoot}{Logical indicating if the results of  estimation  should be printed at each iteration in the  bootstrap restarting  process.}
  \item{pow}{The powers of the pseudo covariates employed by the algorithm.}
  \item{digits}{If specified it means the desidered number of decimal points of the kink estimators to be used during the iterative algorithm.}
  \item{grid}{It measures how close between the two adjacent change points should be merged, default is \code{NULL}.}
  \item{n.boot}{ Positive integer indicating the times of  bootstrap re-sampling in the bootstrap restarting algorithm, default is 20.}
}

\value{
A list with the arguments as components to be used by '\code{mkqr.fit}' and '\code{mkqr.bea}'.
}
\references{
Wei Zhong, Chuang Wan and Wenyang Zhang. (2020) Estimation and inference for multi-kink quantile regression. working paper.
}
\seealso{
\code{\link{mkqr.bea}},\code{\link{mkqr.fit}}
}
\author{
Chuang Wan
}

\examples{
fit.control(K.max=8)
}

