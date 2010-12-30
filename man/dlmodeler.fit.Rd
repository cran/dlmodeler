\name{dlmodeler.fit}
\alias{dlmodeler.fit}
\alias{dlmodeler.fit.MLE}
\alias{dlmodeler.fit.MSE}
\alias{dlmodeler.fit.MAD}
\alias{dlmodeler.fit.sigma}
\title{
Fitting function for a model (MLE, MSE, MAD, sigma)
}
\description{
Fits a DLM by maximum likelihood (MLE), minimum squared errror (MSE),
minimum average deviation (MAD) or minimum standard deviation
(sigma) methods.
}
\usage{
dlmodeler.fit.MLE(yt, build.fun, par,
                  backend = c('KFAS','FKF','dlm'), method = "L-BFGS-B",
                  verbose = FALSE, filter = TRUE, smooth = FALSE,
                  raw.result = FALSE, ...)

dlmodeler.fit.MSE(yt, build.fun, par,
                  ahead, iters = NCOL(yt)-ahead, step = 1, start = 1,
                  backend = c('KFAS','FKF','dlm'), method = "L-BFGS-B",
                  verbose = FALSE, filter = TRUE, smooth = FALSE,
                  raw.result=FALSE, ...)

dlmodeler.fit.MAD(yt, build.fun, par,
                  ahead, iters = NCOL(yt)-ahead, step = 1, start = 1,
                  backend = c('KFAS','FKF','dlm'), method = "L-BFGS-B",
                  verbose = FALSE, filter = TRUE, smooth = FALSE,
                  raw.result=FALSE, ...)

dlmodeler.fit.sigma(yt, build.fun, par,
                  backend = c('KFAS','FKF','dlm'), method = "L-BFGS-B",
                  verbose = FALSE, filter = TRUE, smooth = FALSE,
                  raw.result=FALSE, ...)
}
\arguments{
  \item{yt}{matrix of observed values (one column per time step).}
  \item{build.fun}{function taking parameter vector \code{p} as first argument and returning a DLM.}
  \item{par}{initial value of the parameter vector \code{p}.}
  \item{backend}{an optional argument which specifies the back-end to use for the computations.}
  \item{method}{optimization method passed to function \code{optim}.}
  \item{verbose}{if \code{TRUE}, then write one line per iteration giving the parameter vector \code{p} and the value of the objective function.}
  \item{filter}{if \code{TRUE}, then return the filtered optimal model.}
  \item{smooth}{if \code{TRUE}, the return the smoothed optimal model.}
  \item{raw.result}{if \code{TRUE}, the raw results from the back-end will be stored in \code{raw.result}.}
  \item{ahead}{in case of MSE fitting, the number of predictions to make for each iteration.}
  \item{iters}{in case of MSE fitting, the number of iterations.}
  \item{step}{in case of MSE fitting, the step between iterations.}
  \item{start}{in case of MSE fitting, the index of the first prediction.}
  \item{...}{additional arguments passed to \code{build.fun}.}
}
\details{
\code{dlmodeler.fit.MLE} is designed to find parameter values which maximize the
log-likelihood for the given data. This is called Maximum Likelihood Estimation.

\code{dlmodeler.fit.MSE} is designed to find parameter values which minimize the
average \eqn{n}-step ahead prediction squared error \eqn{(predicted-actual)^2}
for the given data.
This is called Minimum Squared Error fitting. The squared error is averaged over
\code{ahead} prediction steps. Note that having \code{ahead==1} is roughly
equivalent to MLE fitting as long as only the mean is concerned.

\code{dlmodeler.fit.MAD} is designed to find parameter values which minimize the
average \eqn{n}-step ahead prediction absolute error \eqn{|predicted-actual|}
for the given data.
This is called Minimum Average Deviation fitting. The absolute error is averaged over
\code{ahead} prediction steps.

\code{dlmodeler.fit.sigma} is designed to find parameter values which minimize the
one-step ahead prediction variance for the given data.

}
\value{
An object of class \code{dlmodeler.fit} with the following values:
\item{par}{optimal parameter returned by the optimization function \code{optim()}}
\item{message}{message returned by the optimization function \code{optim()}}
\item{convergence}{convergence code returned by the optimization function \code{optim()}}
\item{model}{optimal model found: \code{build.fun(par)}}
\item{logLik}{value of the log-likelihood or NA}
\item{par0}{initial value of \code{par}}
\item{filtered}{optionally, the filtered model: \code{dlmodeler.filter(yt,build.fun(par))}}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Cyrille Szymanski <cnszym@gmail.com>
}
\note{
%%  ~~further notes~~
}

\seealso{
\code{\link{dlmodeler}},
\code{\link{dlmodeler.filter}},
\code{\link{dlmodeler.smooth}},
\code{\link{dlmodeler.forecast}}
}
\examples{
require(dlmodeler)

# analysis from Durbin & Koopman book pp. 32

# load and show the data
y <- matrix(Nile,nrow=1)
plot(y[1,],type='l')

# y(t)   = a(t) + eta(t)
# a(t+1) = a(t) + eps(t)
build.fun <- function(p) {
	varH <- exp(p[1])
	varQ <- exp(p[2])*varH
	dlmodeler.build.polynomial(0,sqrt(varH),sqrt(varQ),name='p32')
}

# fit the model by maximum likelihood estimation
fit <- dlmodeler.fit.MLE(y, build.fun, c(0,0), verbose=FALSE)

# compare the fitted parameters with those reported by the authors
fit$par[2]        # psi = -2.33
fit$model$Ht[1,1] # H   = 15099
fit$model$Qt[1,1] # Q   = 1469.1

# compute the filtered and smoothed values
f <- dlmodeler.filter(y, fit$mod, smooth=TRUE)

# f.ce represents the filtered one steap ahead observation
# prediction expectations E[y(t) | y(1), y(2), ..., y(t-1)]
f.ce <- dlmodeler.extract(f, fit$model,
                          type="observation", value="mean")

# s.ce represents the smoothed observation expectations
# E[y(t) | y(1), y(2), ..., y(n)]
s.ce <- dlmodeler.extract(f$smooth, fit$model,
                          type="observation", value="mean")

# plot the components
plot(y[1,],type='l')
lines(f.ce$p32[1,],col='light blue',lty=2)
lines(s.ce$p32[1,],col='dark blue')
}
\keyword{ dlm }
\keyword{ fit }
\keyword{ MLE }
\keyword{ MSE }
\keyword{ MAD }
\keyword{ sigma }
