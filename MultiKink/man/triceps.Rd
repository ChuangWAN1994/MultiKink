\name{triceps}
\alias{triceps}
\docType{data}
\title{Triceps skinfold thichness dataset}
\description{
The data are derived from an anthropometric study of 892 females under 50  years in three Gambian villages in West Africa.
}
\usage{data(triceps)}
\format{
  A data frame with 892 observations on the following 3 variables.
  \describe{
    \item{\code{age}}{Age of respondents.}
    \item{\code{lntriceps}}{Log of the triceps skinfold thickness.}
    \item{\code{triceps}}{Triceps skinfold thichness.}
  }
}

\source{
Cole T.J., Green P.J. (1992). Smoothing reference centile curves: the LMS method and penalized likelihood. Statistics in medicine, 11(10): 1305-1319.

Perperoglou A, Sauerbrei W, Abrahamowicz M, et al (2019). A review of spline function procedures in R. BMC medical research methodology,  19(1): 46-52.
}
\references{
Cole T.J., Green P.J. (1992). Smoothing reference centile curves: the LMS method and penalized likelihood. Statistics in medicine, 11(10): 1305-1319.

Perperoglou A, Sauerbrei W, Abrahamowicz M, et al (2019). A review of spline function procedures in R. BMC medical research methodology, 19(1): 46-52.
}
\examples{
data(triceps)
## maybe str(triceps) ...
}
\keyword{datasets}
