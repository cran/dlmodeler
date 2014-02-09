\name{dlmodeler.filter.smooth}
\alias{dlmodeler.filter}
\alias{dlmodeler.smooth}
\title{
Filtering and smoothing for a DLM
}
\description{
Kalman filtering and smoothing for a Dynamic Linear Model, using the specified back-end for the computations.
}
\usage{
dlmodeler.filter(yt, model,
                 backend = c("KFAS", "FKF", "dlm"),
                 smooth = FALSE,
                 raw.result = FALSE,
                 logLik = FALSE, filter = TRUE)

dlmodeler.smooth(filt, raw.result = FALSE)
}
\arguments{
  \item{yt}{matrix of observed values (one column per time step).}
  \item{model}{an object of class \code{dlmodeler}.}
  \item{backend}{an optional argument which specifies the back-end to use for the computations.}
  \item{smooth}{an optional argument which specifies if the back-end shoud also run the smoothing algorithm.}
  \item{raw.result}{if TRUE, the raw results from the back-end will be stored in \code{raw.result}.}
  \item{logLik}{an optional argument which specifies if the back-end shoud compute the log-likelihood.}
  \item{filter}{an optional argument which specifies if the back-end shoud also run the filtering algorithm.}
  \item{filt}{filtered \code{dlmodeler.filtered}, as a result from a call to \code{dlmodeler.filter()}.}
}
\details{
This function will automatically load the adequate back-end package.

Currently, packages \code{KFAS} (used by default), \code{FKF} and \code{dlm} are
supported. Refer to \code{\link{dlmodeler}} for more information.
}
\value{
An object of class \code{dlmodeler.filtered} which contains the following elements:
\item{f}{matrix containing the one step ahead predictions \eqn{E(y(t) | y(1),y(2)...y(t-1))}}
\item{at}{matrix containing the one step ahead predicted state variables \eqn{E(a(t) | y(1),y(2)...y(t-1))}}
\item{Pt}{matrix containing the one step ahead predicted state covariance matrices \eqn{cov(a(t) | y(1),y(2)...y(t-1))}}
\item{logLik}{the value of the log-likelihood for the model}
\item{backend}{a character string indicating which back-end was used for the computation}
\item{raw.result}{the raw result from the back-end, or NA if it wasn't requested}

Or an object of class \code{dlmodeler.smoothed} which contains the following elements:
\item{at}{matrix containing the one step ahead smoothed state variables \eqn{E(a(t) | y(1),y(2)...y(n))}}
\item{Pt}{matrix containing the one step ahead predicted state covariance matrices \eqn{cov(a(t) | y(1),y(2)...y(n))}}
\item{backend}{a character string indicating which back-end was used for the computation}
\item{raw.result}{the raw result from the back-end, or NA if it wasn't requested}
}
\author{
Cyrille Szymanski <cnszym@gmail.com>
}
\note{
Package \code{dlm} does not offer a way to obtain the
log-likelihood and the filtered values at the same time (as of v1.1-2).
The log-likelihood is not computed by default, but this can be done by using the
parameter \code{logLik=TRUE}. The computation of the filtered values can also be
disabled with parameter \code{filter=FALSE} if these values are not needed.

Package \code{FKF} does not implement a smoothing algorithm (as of v0.1.1).
}

\seealso{
\code{\link{dlmodeler}}, 
\code{\link{dlmodeler.forecast}}
}
\examples{
\dontrun{
require(dlmodeler)

# generate some data
N <- 365*5
t <- c(1:N,rep(NA,365))
a <- rnorm(N+365,0,.5)
y <- pi + cos(2*pi*t/365.25) + .25*sin(2*pi*t/365.25*3) +
     exp(1)*a + rnorm(N+365,0,.5)

# build a model for this data
m <- dlmodeler.build.polynomial(0,sigmaH=.5,name='level') +
     dlmodeler.build.dseasonal(7,sigmaH=0,name='week') +
     dlmodeler.build.tseasonal(365.25,3,sigmaH=0,name='year') +
     dlmodeler.build.regression(a,sigmaH=0,name='reg')
m$name <- 'mymodel'

system.time(f <- dlmodeler.filter(y, m, raw.result=TRUE))

# extract all the components
m.state.mean <- dlmodeler.extract(f,m,type="state",
                                  value="mean")
m.state.cov <- dlmodeler.extract(f,m,type="state",
                                 value="covariance")
m.obs.mean <- dlmodeler.extract(f,m,type="observation",
                                value="mean")
m.obs.cov <- dlmodeler.extract(f,m,type="observation",
                               value="covariance")
m.obs.int <- dlmodeler.extract(f,m,type="observation",
                               value="interval",prob=.99)

par(mfrow=c(2,1))

# show the one step ahead forecasts & 99\% prediction intervals
plot(y,xlim=c(N-10,N+30))
lines(m.obs.int$mymodel$upper[1,],col='light grey')
lines(m.obs.int$mymodel$lower[1,],col='light grey')
lines(m.obs.int$mymodel$mean[1,],col=2)

# see to which values the filter has converged:
m.state.mean$level[,N] # should be close to pi
mean(abs(m.state.mean$week[,N])) # should be close to 0
m.state.mean$year[1,N] # should be close to 1
m.state.mean$year[6,N] # should be close to .25
m.state.mean$reg[,N] # should be close to e

# show the filtered level+year components
plot(m.obs.mean$level[1,]+m.obs.mean$year[1,],
		type='l',ylim=c(pi-2,pi+2),col='light green',
		ylab="smoothed & filtered level+year")

system.time(s <- dlmodeler.smooth(f,m))

# show the smoothed level+year components
s.obs.mean <- dlmodeler.extract(s,m,type="observation",
                                value="mean")
lines(s.obs.mean$level[1,]+s.obs.mean$year[1,],type='l',
		ylim=c(pi-2,pi+2),col='dark green')
}
}
\keyword{ dlm }
\keyword{ kalman }
\keyword{ filter }
\keyword{ smoother }
