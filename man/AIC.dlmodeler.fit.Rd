\name{AIC.dlmodeler.fit}
\alias{logLik.dlmodeler.filtered}
\alias{logLik.dlmodeler.fit}
\alias{AIC.dlmodeler.fit}
\title{
Log-likelihood and AIC of a model
}
\description{
Returns the log-likelihood or the AIC for a fitted DLM object.
}
\usage{
\method{logLik}{dlmodeler.filtered}(object, ...)

\method{logLik}{dlmodeler.fit}(object, ...)

\method{AIC}{dlmodeler.fit}(object, ..., k = 2)
}
\arguments{
  \item{object}{fitted DLM as given by a call to one of the 
                \code{dlmodeler.fit()} functions, or filtered
                DLM as given by a call to \code{dlmodeler.filter}.}
  \item{...}{not used.}
  \item{k}{penalty parameter.}
}
\details{
The AIC is computed according to the formula
\eqn{-2*log(likelihood) + k*npar},
where \eqn{npar} represents the number of parameters in the fitted model,
and \eqn{k = 2} for the usual AIC, 
or \eqn{k = log(n)} (\eqn{n} the number of observations)
for the BIC or SBC (Schwarz's Bayesian criterion).
}
\value{
Returns a numeric value with the corresponding log-likelihiid, AIC,
BIC, or ..., depending on the value of \code{k}.
}
\references{
Durbin, and Koopman, Time Series Analysis by State Space Methods, Oxford University Press (2001),
page 152.
}
\author{
Cyrille Szymanski <cnszym@gmail.com>
}

\seealso{
\code{\link{dlmodeler.fit.MLE}}
}
\examples{
## Example TODO
}
\keyword{ dlm }
\keyword{ summary }
