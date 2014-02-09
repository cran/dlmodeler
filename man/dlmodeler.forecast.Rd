\name{dlmodeler.forecast}
\alias{dlmodeler.forecast}
\title{
Forecast function
}
\description{
Simulates forecasting for a DLM, with the specified horizon, step and number of iterations.
}
\usage{
dlmodeler.forecast(yt, model,
                   ahead = 1,
                   iters = 1, step = 1, start = 1,
                   prob = .90,
                   backend = c('KFAS','FKF','dlm'),
                   debug = FALSE)
}
\arguments{
  \item{yt}{matrix of observed values (one column per time-step).}
  \item{model}{an instance of \code{dlmodeler}.}
  \item{ahead}{in case of MSE fitting, the number of predictions to make for each iteration.}
  \item{iters}{in case of MSE fitting, the number of iterations.}
  \item{step}{in case of MSE fitting, the step between iterations.}
  \item{start}{in case of MSE fitting, the index of the first prediction.}
  \item{prob}{probability to use for the computation of prediction intervals.}
  \item{backend}{an optional argument which specifies the back-end to use for the computations.}
  \item{debug}{use slow but more robust code.}
}
\details{
This function simulates forecasts for the
specified serie \code{yt} and \code{model}
by generating \code{iters} forecasts every \code{step} points.
The procedure starts at position \code{start}, and each
iteration, \code{ahead} values are predicted.
}
\value{
A \code{data.frame} with the following variables:
\item{index}{the index of the forecasted value, \code{index==i} means that the \eqn{i}-th element of the serie is forecasted}
\item{distance}{the forecasting distance, \code{distance==k} means that this value is a \eqn{k}-step ahead forecast}
\item{lower}{the lower bound for \code{yhat} computed with probability \code{prob}}
\item{yhat}{the forecasted value for the specified \code{index} and \code{distance}}
\item{upper}{the upper bound for \code{yhat} computed with probability \code{prob}}
\item{y}{the observed serie for the specified \code{yt[index]}}
}
\author{
Cyrille Szymanski <cnszym@gmail.com>
}
\note{
Currently, the function only works for univariate time-series.

Currently the implementation is very slow, but its speed will be increased in
future versions of this package.
}

\seealso{
\code{\link{dlmodeler}}, 
\code{\link{dlmodeler.filter}}, 
\code{\link{dlmodeler.smooth}}, 
\code{\link{dlmodeler.forecast}}
}
\examples{
require(dlmodeler)

# generate some quarterly data
n <- 80
level <- 12
sigma <- .75
season <- c(5,6,8,2)
y <- level + 3*sin((1:n)/10) + rep(season,n/4) + rnorm(n, 0, sigma)
y <- matrix(y,nrow=1)

# fit a stochastic level + quarterly seasonal model to the data by
# maximum likelihood estimation
build.fun <- function(p) {
	sigmaH <- exp(p[1])
	sigmaQ <- exp(p[2])*sigmaH
	mod <- dlmodeler.build.polynomial(0,sigmaH=sigmaH,sigmaQ=sigmaQ) +
         dlmodeler.build.dseasonal(4,sigmaH=0)
	return(mod)
}
fit <- dlmodeler.fit.MLE(y, build.fun, c(0,0))

# generate forecasts for observations 81 to 100
f <- dlmodeler.forecast(y, fit$model, start=80, ahead=20)
plot(y[1,],type='l',xlim=c(60,100),ylim=c(10,30))
lines(f$index,f$yhat,col='dark blue')
lines(f$index,f$lower,col='light blue')
lines(f$index,f$upper,col='light blue')

# simulate forecasts post-ex.
f <- dlmodeler.forecast(y, fit$model, ahead=20, start=20, iters=40)
plot(y[1,],type='p')

# show the one step ahead forecasts
with(f[f$distance==1,], lines(index,yhat,col='dark blue'))
# show the 10 step ahead forecasts
with(f[f$distance==10,], lines(index,yhat,col='blue'))
# show the 20 step ahead forecasts
with(f[f$distance==20,], lines(index,yhat,col='light blue'))

}
\keyword{ dlm }
\keyword{ kalman }
\keyword{ forecast }
