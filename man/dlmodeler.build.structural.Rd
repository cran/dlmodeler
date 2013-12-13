\name{dlmodeler.build.structural}
\alias{dlmodeler.build.structural}
\title{
Build a structural time series model
}
\description{
Builds a DLM for a structural time series, consisting of
a polynomial term (level, trend, ...),
a "dummy seasonal" pattern,
a trigonometric cycle term,
and an observation disturbance.
}
\usage{
dlmodeler.build.structural(pol.order = NULL, dseas.order = NULL,
                           tseas.period = NULL, tseas.order = NULL,
                           sigmaH = NA, pol.sigmaQ = 0,
                           dseas.sigmaQ = 0, tseas.sigmaQ = 0,
                           name = "structural")
}
\arguments{
  \item{pol.order}{order of the polynomial (0=constant, 1=linear, 2=cubic...), or \code{NULL} to ignore polynomial component.}
  \item{dseas.order}{period of the dummy seasonal pattern, or \code{NULL} to ignore dummy seasonal component.}
  \item{tseas.period}{period of the trigonometric seasonal pattern, or \code{NULL} to ignore trigonometric component.}
  \item{tseas.order}{number of harmonics in the trigonometric seasonal pattern, or \code{NULL}. Optional when tseas.period is an integer (a default value is used), mandatory otherwise.}
  \item{sigmaH}{std dev of the observation disturbance (if unknown, set to NA and use dlmodeler.fit to estimate it).}
  \item{pol.sigmaQ}{std dev of the polynomial state disturbances (if unknown, set to NA and use dlmodeler.fit to estimate it).}
  \item{dseas.sigmaQ}{std dev of the dummy seasonal state disturbances (if unknown, set to NA and use dlmodeler.fit to estimate it).}
  \item{tseas.sigmaQ}{std dev of the trigonometric seasonal state disturbances (if unknown, set to NA and use dlmodeler.fit to estimate it).}
  \item{name}{an optional name to be given to the resulting DLM.}
}
\details{
The initial value \code{P0inf} is parametered to use exact diffuse initialisation
(if supported by the back-end). 
}
\value{
An object of class \code{dlmodeler} representing the structural model.
This object can have the following components:
\item{level}{component representing the level (when \code{pol.order} = 0)}
\item{level+trend}{component representing the level+trend (when \code{pol.order} = 1)}
\item{polynomial}{component representing the level, trend, ... (when \code{pol.order} > 1)}
\item{seasonal}{component representing the dummy seasonal pattern}
\item{trigonometric}{component representing the trigonometric seasonal pattern}
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
\code{\link{dlmodeler.build.dseasonal}}, 
\code{\link{dlmodeler.build.tseasonal}}, 
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

# deterministic level + quarterly seasonal
mod <- dlmodeler.build.structural(pol.order=0, dseas.order=4,
                                  sigmaH=sigma)
f <- dlmodeler.filter(y, mod)

# show the one step ahead forecasts
par(mfrow=c(2,1))
plot(y,type='l')
lines(f$f[1,],col='light blue')

# show the filtered level and seasonal components
c <- dlmodeler.extract(f,mod,type="state")
lines(c$level[1,],col='blue')
plot(c$seasonal[1,],type='l',col='dark green')
}
}
\keyword{ dlm }
\keyword{ structural }
