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
                           tseas.order = NULL, tseas.period = NULL,
                           sigmaH = NA, pol.sigmaQ = 0,
                           dseas.sigmaQ = 0, tseas.sigmaQ = 0,
                           name = "structural")
}
\arguments{
  \item{pol.order}{order of the polynomial (0=constant, 1=linear, ...), or \code{NULL}.}
  \item{dseas.order}{period of the dummy seasonal pattern, or \code{NULL}.}
  \item{tseas.order}{number of harmonics in the trigonometric seasonal pattern, or \code{NULL}.}
  \item{tseas.period}{period of the trigonometric seasonal pattern, or \code{NULL}.}
  \item{sigmaH}{std dev of the observation disturbance.}
  \item{pol.sigmaQ}{std dev of the polynomial state disturbances.}
  \item{dseas.sigmaQ}{std dev of the dummy seasonal state disturbances.}
  \item{tseas.sigmaQ}{std dev of the trigonometric seasonal state disturbances.}
  \item{name}{an optional name to be given to the resulting DLM.}
}
\details{
The initial value \code{P0inf} is parametered to use exact diffuse initialisation
(if supported by the back-end). 

Remember: all functions of this package assume by default a zero state
disturbance (deterministic states). By looking at the prototype
for this function, you can assume sigmaH refers to the observation and sigmaQ
refers to the state. This rule of thumb is useful to remember how this function
should be used.
}
\value{
An object of class \code{dlmodeler} representing the structural model.
This object can have the following components:
\item{level}{component representing the level (when \code{pol.order} = 0)}
\item{level+trend}{component representing the level+trend (when \code{pol.order} = 1)}
\item{polynomial}{component representing the level, trend, ... (when \code{pol.order} > 1)}
\item{seasonal}{component representing the dummy seasonal pattern}
\item{cycle}{component representing the trigonometric seasonal cycle}
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
\keyword{ dlm }
\keyword{ structural }
