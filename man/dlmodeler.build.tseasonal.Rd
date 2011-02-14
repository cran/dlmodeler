\name{dlmodeler.build.tseasonal}
\alias{dlmodeler.build.tseasonal}
\title{
Build a trigonometric seasonal model
}
\description{
Builds an univariate trigonometric seasonal DLM of the specified order.
}
\usage{
dlmodeler.build.tseasonal(per, ord = NULL,
                          sigmaH = NA, sigmaQ = 0,
                          name = "tseasonal")
}
\arguments{
  \item{per}{period of the seasonal pattern.}
  \item{ord}{optional order (number of harmonics) of the seasonal pattern.}
  \item{sigmaH}{std dev of the observation disturbance.}
  \item{sigmaQ}{std dev of the state disturbance.}
  \item{name}{an optional name to be given to the resulting DLM.}
}
\details{
The trigonometric decomposition has the form
\eqn{a[1]cos(2pi/per) + a[2]sin(2pi/per) + a[3]cos(2pi/per*2) + a[4]sin(2pi/per*2) ... + a[2*ord-1]cos(2pi/per*ord) + a[2*ord]sin(2pi/per*ord)}.

If \code{ord} is not specified, the order is selected such that there are
\code{per}-1 coefficients in the decomposition. In this case, \code{per} must
be an integer value.

The initial value \code{P0inf} is parametered to use exact diffuse initialisation
(if supported by the back-end). 

Remember: all functions of this package assume by default a zero state
disturbance (deterministic states). By looking at the prototype
for this function, you can assume sigmaH refers to the observation and sigmaQ
refers to the state. This rule of thumb is useful to remember how this function
should be used.
}
\value{
An object of class \code{dlmodeler} representing the trigonometric seasonal model.
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
\code{\link{dlmodeler.build.structural}}, 
\code{\link{dlmodeler.build.arima}}, 
\code{\link{dlmodeler.build.regression}}
}
\examples{
require(dlmodeler)

# generate some data
N <- 365*5
t <- c(1:N,rep(NA,365))
a <- rnorm(N+365,0,.5)
y <- pi + cos(2*pi*t/365.25) + .25*sin(2*pi*t/365.25*3) +
     exp(1)*a + rnorm(N+365,0,.5)

# build a model for this data
m1 <- dlmodeler.build.polynomial(0,sigmaH=.5,name='level')
m2 <- dlmodeler.build.dseasonal(7,sigmaH=0,name='week')
m3 <- dlmodeler.build.tseasonal(365.25,3,sigmaH=0,name='year')
m4 <- dlmodeler.build.regression(a,sigmaH=0,name='reg')
m <- dlmodeler.add(m1,dlmodeler.add(m2,dlmodeler.add(m3,m4)),
                   name='mymodel')

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
\keyword{ dlm }
\keyword{ seasonal }
\keyword{ trigonometric }
