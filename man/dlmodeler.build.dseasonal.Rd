\name{dlmodeler.build.dseasonal}
\alias{dlmodeler.build.dseasonal}
\title{
Build a "dummy seasonal" model
}
\description{
Builds an univariate "dummy seasonal" DLM of the specified order.
}
\usage{
dlmodeler.build.dseasonal(ord, sigmaH = NA, sigmaQ = 0,
                          name = "dseasonal")
}
\arguments{
  \item{ord}{period of the seasonal pattern.}
  \item{sigmaH}{std dev of the observation disturbance (if unknown, set to NA and use dlmodeler.fit to estimate it).}
  \item{sigmaQ}{std dev of the state disturbance (if unknown, set to NA and use dlmodeler.fit to estimate it).}
  \item{name}{an optional name to be given to the resulting DLM.}
}
\details{
The seasonal pattern is represented by \code{ord} seasonal indices
\eqn{a[1], a[2], ..., a[ord]}.
The indices are constrained such that their sum equals 0, with
\eqn{a[ord] = -a[1] - a[2] - a[3] ... -a[ord-1]}.
This only requires \code{ord}-1 state variables.

The initial value \code{P0inf} is parametered to use exact diffuse initialisation
(if supported by the back-end).
}
\value{
An object of class \code{dlmodeler} representing the dummy seasonal model.
}
\references{
Durbin, and Koopman, Time Series Analysis by State Space Methods, Oxford University Press (2001),
pages 38-45.
}
\author{
Cyrille Szymanski <cnszym@gmail.com>
}
\note{
State representations are not unique, so other forms could be used to achieve the same goals.
}

\seealso{
\code{\link{dlmodeler}}, 
\code{\link{dlmodeler.build}}, 
\code{\link{dlmodeler.build.polynomial}}, 
\code{\link{dlmodeler.build.tseasonal}}, 
\code{\link{dlmodeler.build.structural}}, 
\code{\link{dlmodeler.build.arima}}, 
\code{\link{dlmodeler.build.regression}}
}
\examples{
\dontrun{
require(dlmodeler)

# generate some quarterly data
n <- 80
level <- 12
sigma <- .75
season <- c(5,6,8,2)
y <- level + rep(season,n/4) + rnorm(n, mean=0, sd=sigma)

# deterministic level + quarterly seasonal + disturbance
mod1 <- dlmodeler.build.polynomial(0,sigmaH=sigma)
mod2 <- dlmodeler.build.dseasonal(4,sigmaH=0)
mod <- dlmodeler.add(mod1, mod2)
f <- dlmodeler.filter(y, mod)

# show the one step ahead forecasts
plot(y,type='l')
lines(f$f[1,],col='light blue')
}
}
\keyword{ dlm }
\keyword{ seasonal }
\keyword{ dummy }
