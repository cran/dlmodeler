\name{AIC.dlmodeler.fit}
\alias{AIC.dlmodeler.fit}
\title{
AIC of a model
}
\description{
Returns the AIC for a fitted DLM object.
}
\usage{
\method{AIC}{dlmodeler.fit}(object, ..., k = 2)
}
\arguments{
  \item{object}{fitted DLM as given by a call to one of the 
                \code{dlmodeler.fit()} functions.}
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
Returns a numeric value with the corresponding AIC, BIC, or ...,
depending on the value of \code{k}.
}
\references{
Durbin, and Koopman, Time Series Analysis by State Space Methods, Oxford University Press (2001),
pp. 152.
}
\author{
Cyrille Szymanski <cnszym@gmail.com>
}
\note{
%%  ~~further notes~~
}

\seealso{
\code{\link{dlmodeler.fit.MLE}}
}
\examples{
## Example TODO
}
\keyword{ dlm }
\keyword{ summary }
