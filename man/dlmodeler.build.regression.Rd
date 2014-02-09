\name{dlmodeler.build.regression}
\alias{dlmodeler.build.regression}
\alias{dlmodeler.regression}
\title{
Build a regression model
}
\description{
Builds an univariate (multi-linear) regression DLM.
}
\usage{
dlmodeler.regression(covariates, sigmaH = NA, sigmaQ = 0,
                     intercept = FALSE, name = "regression")

dlmodeler.build.regression(covariates, sigmaH = NA, sigmaQ = 0,
                           intercept = FALSE, name = "regression")
}
\arguments{
  \item{covariates}{covariate matrix (one row per covariate).}
  \item{sigmaH}{std dev of the observation disturbance (if unknown, set to NA and use dlmodeler.fit to estimate it). Default = NA.}
  \item{sigmaQ}{std dev of the state disturbance (if unknown, set to NA and use dlmodeler.fit to estimate it). Default = 0.}
  \item{intercept}{should an intercept be added to the model?}
  \item{name}{an optional name to be given to the resulting DLM.}
}
\details{
The regression term is of the form
\eqn{a[1] x_1(t) + a[2]x_2(t) ... + a[k]x_k(t)}, where
\eqn{x_k} is the \eqn{k}-th covariate.

The initial value \code{P0inf} is parametered to use exact diffuse initialisation
(if supported by the back-end).
}
\value{
An object of class \code{dlmodeler} representing the regression model.
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
\code{\link{dlmodeler.build.structural}}, 
\code{\link{dlmodeler.build.arima}}
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
     dlmodeler.build.dseasonal(7,sigmaH=0,name='week')
     dlmodeler.build.tseasonal(365.25,3,sigmaH=0,name='year')
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
\keyword{ regression }
