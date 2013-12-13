\name{dlmodeler.build.arima}
\alias{dlmodeler.build.arima}
\title{
Build an ARIMA model
}
\description{
Builds an univariate ARIMA DLM of the specified order and coefficients.
}
\usage{
dlmodeler.build.arima(ar=c(), ma=c(), d=0,
                      sigmaH = NA, sigmaQ = 0,
                      name = "arima")
}
\arguments{
  \item{ar}{vector of autoregressive coefficients c(ar1, ar2, ar3...).}
  \item{ma}{vector of moving average coefficients c(ma1, ma2, ma3...).}
  \item{d}{order of differenciation.}
  \item{sigmaH}{std dev of the observation disturbance (if unknown, set to NA and use dlmodeler.fit to estimate it).}
  \item{sigmaQ}{std dev of the state disturbances (if unknown, set to NA and use dlmodeler.fit to estimate it).}
  \item{name}{an optional name to be given to the resulting DLM.}
}
\details{
The autoregressive terms of the model are \eqn{ar[1] + ar[2]L + ... ar[p]L^p}
where \eqn{L} is the lag operator.

The moving average terms of the model are \eqn{1 + ma[1]L + ... ma[q]L^q}
where \eqn{L} is the lag operator.

The initial value \code{P0inf} is parametered to use exact diffuse initialisation
(if supported by the back-end).
}
\value{
An object of class \code{dlmodeler} representing the ARIMA model.
}
\references{
Durbin, and Koopman, Time Series Analysis by State Space Methods, Oxford University Press (2001),
pages 46-48.
}
\author{
Cyrille Szymanski <cnszym@gmail.com>
}
\note{
State representations are not unique, so other forms could be used to achieve the same goals.

Currently, only ARMA models (d=0) are implemented.
}

\seealso{
\code{\link{dlmodeler}}, 
\code{\link{dlmodeler.build}}, 
\code{\link{dlmodeler.build.polynomial}}, 
\code{\link{dlmodeler.build.dseasonal}}, 
\code{\link{dlmodeler.build.tseasonal}}, 
\code{\link{dlmodeler.build.structural}}, 
\code{\link{dlmodeler.build.regression}}
}
\examples{
# Example TODO
}
\keyword{ dlm }
\keyword{ arima }
