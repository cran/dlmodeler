\name{dlmodeler.extract}
\alias{dlmodeler.extract}
\title{
Extract the mean, covariance and prediction intervals for states and observations
}
\description{
Extracts the mean (expectation), the variance-covariance matrix,
and the prediction intervals for the states and observations of a
filtered or smoothed DLM component.
}
\usage{
dlmodeler.extract(fs, model, compnames=NULL,
                  type=c("observation","state"),
                  value=c("mean","covariance","interval"),
                  prob=.90)
}
\arguments{
  \item{fs}{filtered or smoothed \code{dlmodeler}, as a result from a call
  			to \code{dlmodeler.filter()} or \code{dlmodeler.smooth()}.}
  \item{model}{object of class \code{dlmodeler} which was used
               for filtering or smoothing.}
  \item{compnames}{an optional list of components to extract.}
  \item{type}{an optional string indicating the type to extract:
              observation (output, by default) or state.}
  \item{value}{an optional string indicating the value to extract:
               mean (expectation, by default), covariance matrix,
               or prediction intervals.}
  \item{prob}{an optional probability (default = 90\%) for the
              computation of prediction intervals.}
}
\details{
A component is a named portion of the state vector matrix which can be
extracted with this function.
Components are automatically created when DLMs are added together
which makes it easier to decompose it later into its building blocks
(for example: level+trend+seasonal+cycle).

Let us assume model named \code{m} is constructed by adding models
named \code{m1} and \code{m2}. Typically, \code{m} will be constructed
with two components named \code{m1} and \code{m1}, which can be
extracted by this function.
}
\value{
When this function is used with a filtered \code{dlmodeler}, it returns
the means and covariances of the one-step ahead forecasts for the components:
\itemize{
\item \code{Zt \%*\% at} \eqn{= E(y(t) | y(1),y(2)...y(t-1))} for observation means, in the form of a \eqn{(d,n)} matrix.
\item \code{at} \eqn{= E(alpha(t) | y(1),y(2)...y(t-1))} for state means, in the form of a \eqn{(m,n)} matrix.
\item \code{Zt \%*\% Pt \%*\% t(Zt) + Ht} \eqn{= cov(y(t) | y(1),y(2)...y(t-1))} for observation covariances, in the form of a \eqn{(d,d,n)} array.
\item \code{Pt} \eqn{= cov(alpha(t) | y(1),y(2)...y(t-1))} for state covariances, in the form of a \eqn{(m,m,n)} array.
}

When this function is used with a smoothed \code{dlmodeler}, it returns
the means and covariances of the smoothed components:
\itemize{
\item \code{Zt \%*\% at} \eqn{= E(y(t) | y(1),y(2)...y(N))} for observation means, in the form of a \eqn{(d,n)} matrix.
\item \code{at} \eqn{= E(alpha(t) | y(1),y(2)...y(N))} for state means, in the form of a \eqn{(m,n)} matrix.
\item \code{Zt \%*\% Pt \%*\% t(Zt) + Ht} \eqn{= cov(y(t) | y(1),y(2)...y(N))} for observation covariances, in the form of a \eqn{(d,d,n)} array.
\item \code{Pt} \eqn{= cov(alpha(t) | y(1),y(2)...y(N))} for state covariances, in the form of a \eqn{(m,m,n)} array.
}

When the value \code{interval} is requested, this function returns a list for each component containing:
\itemize{
\item \code{mean} = the mean (expectaton) for the filtered or smoothed state or observation variable.
\item \code{lower} = lower bound of the prediction interval computed as \code{mean-k*sd},
	\code{k=-qnorm((1+prob)/2)}.
\item \code{upper} = upper bound of the prediction interval computed as \code{mean+k*sd},
	\code{k=-qnorm((1+prob)/2)}.
}
}
\author{
Cyrille Szymanski <cnszym@gmail.com>
}

\seealso{
\code{\link{dlmodeler}}, 
\code{\link{dlmodeler.filter}}, 
\code{\link{dlmodeler.smooth}}
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
}
\keyword{ dlm }
\keyword{ components }
