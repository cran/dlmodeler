\name{dlmodeler-package}
\alias{dlmodeler}
\docType{package}
\title{
Generalized Dynamic Linear Modeler
}
\description{
Package \code{dlmodeler} is a set of user friendly functions
to simplify the state-space modelling, fitting, analysis and
forecasting of Generalized Dynamic Linear Models (DLMs).

It includes functions to name and extract components
of a DLM, and provides a unified interface
compatible with other state-space packages including: \code{KFAS},
\code{FKF} and \code{dlm}.
}
\details{
Generalized Dynamic Linear Models are a powerful approach
to time-series modelling, analysis and forecasting. This framework
is closely related to the families of ARMA models, exponential
smoothing, and structural time-series (also known as unobserved
component models) which it includes as special cases.

The origin of DLM time-series analysis has its roots in the
world of engineering. In order to control dynamic physical systems,
unknown quantities such as velocity and position (the state of the
system) need to be estimated from noisy measurements such as readings
from various sensors (the observations). The state of the
system evolves from one state to another according to a defined
transition equation (possibly including random perturbations and
intervention effects), and the observations are derived from the
state values by a defined observation equation (also possibly
including random disturbances and intervention effects).

The challenge is to obtain the best estimate of the unknown state
considering the set of available observations at a given point in
time. During the 1960s, the Kalman filtering and smoothing algorithm
was developed and popularized to efficiently and optimally solve this
etimation problem.
The technique is based on an iterative procedure in which state values
are successively predicted given the knowledge of the past
observations, and then updated upon the reception of the next
observation. Because of the predict-and-update nature of Kalman
filtering, it can also be interpreted under a Bayesian perspective.

The theory developed for the control of dynamic systems has a direct
application to the general analysis of time-series. By having a good
estimate of the current state and dynamics of the system, it is
possible to derive assumptions about their evolution and subsequent
values; and therefore to obtain a forecast for the future
observations.

Dynamic Linear Models are a special case of general state-space
models where the state and the observation equations are linear, and
the distributions follow a normal law. They are also referred to as
gaussian linear state-space models. Generalized DLMs relax the
assumption of normality by allowing the distribution to be any of the
exponential family of functions (which includes the Bernoulli,
binomial and Poisson distributions, useful in particular for count
data).

The distinguishing aspect of this package is that it provides
functions for naming and extracting components of DLMs (see below),
and a unified interface compatible with other state-space packages for
filtering and smoothing:
\itemize{
\item package \code{KFAS}: implements exact diffuse initialization and
      supports filtering, smoothing and likelihood computation for the
      exponential family state-space models (as of v0.6.0)
\item package \code{FKF}: very fast and memory efficient but has no
      smoothing algorithm (as of v0.1.1)
\item package \code{dlm}: good general purpose package with many
      helper functions available and some support for Bayesian
      analysis (as of v1.1-2)
}

}
\section{State-space form and notations}{
The state space model is represented as follows:

initial state:        \eqn{alpha(0) \sim N(a(0), P(0))}

observation equation: \eqn{y(t)       = Z(t)alpha(t) + eta(t)}

transition equation:  \eqn{alpha(t+1) = T(t)alpha(t) + R(t)eps(t)}

observation disturbance: \eqn{eta(t) \sim N(0, H(t))}

state disturbance:       \eqn{eps(t) \sim N(0, Q(t))}

state mean and covariance matrix: \eqn{a(t) = E[alpha(t)]}
and \eqn{P(t) = cov(alpha(t))}

With:
\itemize{
\item \eqn{n} = number of time-steps \eqn{(t=1..n)}
\item \eqn{m} = dimension of state vector
\item \eqn{alpha(t)} = state vector \eqn{(m,1)}
\item \eqn{a(0)} = initial state vector \eqn{(m,1)}
\item \eqn{P(0)} = initial state covariance matrix \eqn{(m,m)}
\item \eqn{Pinf(0)} = diffuse part of \eqn{P(0)} matrix \eqn{(m,m)}
\item \eqn{d} = dimension of observation vector
\item \eqn{y(t)} = observation vector \eqn{(d,1)}
\item \eqn{Z(t)} = observation design matrix \eqn{(d,m)} if constant, or \eqn{(d,m,n)}
\item \eqn{eta(t)} = observation disturbance vector \eqn{(d,1)}
\item \eqn{H(t)} = observation disturbance covariance matrix \eqn{(d,d)} if constant, or \eqn{(d,d,n)}
\item \eqn{T(t)} = state transition matrix \eqn{(m,m)} if constant, or \eqn{(m,m,n)}
\item \eqn{r} = dimension of the state disturbance covariance matrix
\item \eqn{eps(t)} = state disturbance vector \eqn{(r,1)}
\item \eqn{R(t)} = state disturbance selection matrix \eqn{(m,r)} if constant, or \eqn{(m,r,n)}
\item \eqn{Q(t)} = state disturbance covariance matrix \eqn{(r,r)} if constant, or \eqn{(r,r,n)}
}
}
\section{Components}{
DLMs are constructed by combining several terms together. The model
consisiting of \eqn{level+trend+seasonal+cycle} is a good example of
how individual elements can be added together to form a more complete
model. A typical analysis will consider the model as a whole, but also
look into the values of indivdual terms, for example the variations in
the \eqn{level}, and the evolution of the shape of the \eqn{seasonal}
terms.

This package introduces a notion called the "component" to facilitate
the analysis of the DLMs. Mathematically speaking, a component is a
named subset of state variables.

Components can be created when the model is built, which makes it
easier to access and analyze their values afterwards.
}
\section{Notes}{
Work is in progress to include package \code{sspir} as back-end and to
provide support for generalized models.

Ideas for the future include: using formulas to specify DLMs and
faster C code whenever possible.
}
\author{
Cyrille Szymanski <cnszym@gmail.com>
}
\section{Maintainer}{
Cyrille Szymanski <cnszym@gmail.com>
}
\references{
Harvey, A.C. Forecasting, Structural Time Series Models and the Kalman Filter. Cambridge University Press (1989).

Durbin, J. and Koopman, S. J. Time Series Analysis by State Space Methods. Oxford University Press (2001). http://www.ssfpack.com/dkbook/

Commandeur, and Koopman, An Introduction to State Space Time Series Analysis, Oxford University Press (2007).

Petris, Petrone, and Campagnoli, Dynamic Linear Models with R, Springer (2009).
}
\keyword{ dlm }
\keyword{ forecast }
\keyword{ kalman }
\keyword{ filter }
\keyword{ smoother }
\seealso{
Other R packages and functions of interest (in alphabetical order):
\itemize{
\item Package \code{dlm}
\item Package \code{dse}
\item Package \code{FKF}
\item Package \code{forecast}
\item HoltWinters from \{stats\}
\item Package \code{KFAS}
\item Package \code{MARSS}
\item Package \code{sspir}
\item StructTS from \{stats\}
}

Other software of interest (in alphabetical order):
\itemize{
\item SAS PROC UCM \code{http://www.sas.com/products/ets}
\item SsfPack \code{http://www.ssfpack.com}
\item STAMP \code{http://www.stamp-software.com}
}
}
\examples{
require(dlmodeler)

# This section illustrates most of the possibilities offered by the
# package by reproducing famous examples of state-space time-series
# analysis found in the litterature.

##############################################
# analysis from Durbin & Koopman book pp. 32 #
# random walk                                #
##############################################

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


###############################################
# analysis from Durbin & Koopman book pp. 163 #
# random walk + stochastic seasonal           #
###############################################

# load and show the data
y <- matrix(log(Seatbelts[,'drivers']),nrow=1)
plot(y[1,],type='l')

# y(t)    = a(t) + s1(t) + eta(t)
# a(t+1)  = a(t) + eps_L(t)
# s1(t+1) = -s2(t) - s3(t) - ... - s12(t) + eps_S(t)
# s2(t+1) = s1(t)
# s3(t+1) = s2(t), etc.
build.fun <- function(p) {
	sigmaH <- exp(p[1])
	sigmaQL <- exp(p[2])*sigmaH
	sigmaQS <- exp(p[3])*sigmaH
	dlmodeler.build.structural(
			pol.order=0,
			pol.sigmaQ=sigmaQL,
			dseas.order=12,
			dseas.sigmaQ=sigmaQS,
			sigmaH=sigmaH,
			name='p163')
}

# fit the model by maximum likelihood estimation
fit <- dlmodeler.fit.MLE(y, build.fun, c(0,0,0), verbose=FALSE)

# compare the fitted parameters with those reported by the authors
fit$model$Ht[1,1] # H  = 0.0034160
fit$model$Qt[1,1] # Q1 = 0.00093585
fit$model$Qt[2,2] # Q2 = 5.0109e-007

# compute the filtered and smoothed values
f <- dlmodeler.filter(y, fit$model, smooth=TRUE)

# f.ce represents the filtered one steap ahead observation
# prediction expectations E[y(t) | y(1), y(2), ..., y(t-1)]
f.ce <- dlmodeler.extract(f, fit$model,
                          type="observation", value="mean")

# s.ce represents the smoothed observation expectations
# E[y(t) | y(1), y(2), ..., y(n)]
s.ce <- dlmodeler.extract(f$smooth, fit$model,
                          type="observation", value="mean")

# plot the components
plot(y[1,])
lines(f.ce$level[1,],col='light blue',lty=2)
lines(s.ce$level[1,],col='dark blue')

# note that the smoothed seasonal component appears to be constant
# throughout the serie, this is due to Qt[2,2]=sigmaQS being close to
# zero. Durbin & Koopman treat the seasonal component as deterministic
# in the remainder of their models.
plot(y[1,]-s.ce$level[1,],ylim=c(-.5,.5))
lines(f.ce$seasonal[1,],type='l',col='light green',lty=2)
lines(s.ce$seasonal[1,],type='l',col='dark green')


#########################################################
# analysis from Durbin & Koopman book pp. 166           #
# random walk + seasonal + seat belt law + petrol price #
#########################################################

# load and show the data
y <- matrix(log(Seatbelts[,'drivers']),nrow=1)
law <- matrix(Seatbelts[,'law'],nrow=1)
petrolprice <- matrix(log(Seatbelts[,'PetrolPrice']),nrow=1)
par(mfrow=c(3,1))
plot(y[1,],type='l')
plot(petrolprice[1,],type='l')
plot(law[1,],type='l')

# y(t)    = a(t) + s1(t) + lambda*law + mu*petrolprice + eta(t)
# a(t+1)  = a(t) + eps(t)
# s1(t+1) = -s2(t) - s3(t) - ... - s12(t)
# s2(t+1) = s1(t)
# s3(t+1) = s2(t), etc.
build.fun <- function(p) {
	sigmaH <- exp(p[1])
	sigmaQ <- exp(p[2])*sigmaH
	m1 <- dlmodeler.build.structural(
			pol.order=0,
			dseas.order=12,
			sigmaH=sigmaH,
			pol.sigmaQ=sigmaQ)
	m2 <- dlmodeler.build.regression(
			law,
			sigmaH=0,
			name='law')
	m3 <- dlmodeler.build.regression(
			petrolprice,
			sigmaH=0,
			name='petrolprice')
	dlmodeler.add(dlmodeler.add(m1,m2),m3,name='p166')
}

# fit the model by maximum likelihood estimation
fit <- dlmodeler.fit.MLE(y, build.fun, c(0,0), verbose=FALSE)

# compute the filtered and smoothed values
f <- dlmodeler.filter(y, fit$model, smooth=TRUE)

# E[y(t) | y(1), y(2), ..., y(t-1)]
f.ce <- dlmodeler.extract(f, fit$model,
                          type="observation", value="mean")
# E[y(t) | y(1), y(2), ..., y(n)]
s.ce <- dlmodeler.extract(f$smooth, fit$model,
                          type="observation", value="mean")
# E[a(t) | y(1), y(2), ..., y(t-1)]
fa.ce <- dlmodeler.extract(f, fit$model,
                           type="state", value="mean")
# E[a(t) | y(1), y(2), ..., y(n)]
sa.ce <- dlmodeler.extract(f$smooth, fit$model,
                           type="state", value="mean")

# see to which values the model has converged and
# compare them with those reported by the authors
fa.ce$law[1,193]         # law         = -0.23773
fa.ce$petrolprice[1,193] # petrolprice = -0.29140

# plot the smoothed de-seasonalized serie
par(mfrow=c(1,1))
plot(y[1,])
lines(s.ce$level[1,]+s.ce$law[1,]+s.ce$petrolprice[1,],
      col='dark blue')

# show the AIC of the model
AIC(fit)


###################################
# testing other fitting functions #
###################################

# load and show the data
y <- matrix(Nile,nrow=1)
plot(y[1,],type='l')

# random walk
# y(t)   = a(t) + eta(t)
# a(t+1) = a(t) + eps(t)
build.fun <- function(p) {
	if( length(p)==1 ) p <- c(0,p[1])
	sigmaH <- exp(p[1])
	sigmaQ <- exp(p[2])*sigmaH
	dlmodeler.build.polynomial(0,sigmaH,sigmaQ,name='p32')
}

# fit the model by maximum likelihood estimation and compute the
# 1-step ahead MSE
fit.mle <- dlmodeler.fit.MLE(y, build.fun, c(0,0))
mean((fit.mle$filtered$f[1,1:100]-y[1,1:100])^2)

# fit the model by minimizing the 1-step ahead MSE
fit.mse1 <- dlmodeler.fit.MSE(y, build.fun, ahead=1, c(0))
mean((fit.mse1$filtered$f[1,1:100]-y[1,1:100])^2)

# fit the model by minimizing the 3-step ahead MSE
fit.mse3 <- dlmodeler.fit.MSE(y, build.fun, ahead=3, c(0))
mean((fit.mse3$filtered$f[1,1:100]-y[1,1:100])^2)

# fit the model by minimizing the std deviation of the 1-step
# ahead predictions
fit.sigma <- dlmodeler.fit.sigma(y, build.fun, c(0))
mean((fit.sigma$filtered$f[1,1:100]-y[1,1:100])^2)

# compare the 1-step ahead forecasts for these models
# as can be expected, the MLE and MSE1 models roughly
# have the same means
plot(y[1,],type='l')
lines(fit.mle$filtered$f[1,],col='dark blue')
lines(fit.mse1$filtered$f[1,],col='dark green')
lines(fit.mse3$filtered$f[1,],col='dark red')
lines(fit.sigma$filtered$f[1,],col='orange')


#################################################
# looking at variances and prediction intervals #
#################################################

# load and show the data
y <- matrix(log(Seatbelts[,'drivers']),nrow=1)
plot(y[1,],type='l')

# model with level + seasonal
build.fun <- function(p) {
	sigmaH <- exp(p[1])
	sigmaQL <- exp(p[2])*sigmaH
	sigmaQS <- exp(p[3])*sigmaH
	dlmodeler.build.structural(
			pol.order=0,
			pol.sigmaQ=sigmaQL,
			dseas.order=12,
			dseas.sigmaQ=sigmaQS,
			sigmaH=sigmaH,
			name='p163')
}

# fit the model by maximum likelihood estimation, filter & smooth
fit <- dlmodeler.fit.MLE(y, build.fun, c(0,0,0), verbose=FALSE)
fs <- dlmodeler.filter(y, fit$model, smooth=TRUE)

# value we will be using to compute 90\% prediction intervals
prob <- 0.90

# true output mean + prediction interval
output.intervals <- dlmodeler.extract(fs,fit$model,
                    type="observation",value="interval",prob=prob)
plot(y[1,],xlim=c(100,150))
lines(output.intervals$p163$mean[1,],col='dark green')
lines(output.intervals$p163$lower[1,],col='dark grey')
lines(output.intervals$p163$upper[1,],col='dark grey')

# true state level mean + prediction interval
state.intervals <- dlmodeler.extract(fs, fit$model,type="state",
                                     value="interval",prob=prob)
plot(y[1,])
lines(state.intervals$level$mean[1,],col='dark green')
lines(state.intervals$level$lower[1,],col='dark grey')
lines(state.intervals$level$upper[1,],col='dark grey')

# true state seasonal mean + prediction interval
plot(state.intervals$seasonal$mean[1,],
     ylim=c(-.4,.4),xlim=c(100,150),
     type='l',col='dark green')
lines(state.intervals$seasonal$lower[1,],col='light grey')
lines(state.intervals$seasonal$upper[1,],col='light grey')
}
