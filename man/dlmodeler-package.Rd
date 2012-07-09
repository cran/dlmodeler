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
compatible with other state-space packages including: \code{dlm},
\code{KFAS}, \code{FKF} and \code{sspir}.
}
\section{Introduction}{
Generalized Dynamic Linear Models are a powerful approach
to time-series modelling, analysis and forecasting. This framework
is closely related to the families of regression models, ARIMA models,
exponential smoothing, and structural time-series (also known as
unobserved component models, UCM).

The origin of DLM time-series analysis has its roots in the
world of engineering. In order to control dynamic physical systems,
unknown quantities such as velocity and position (the state of the
system) need to be estimated from noisy measurements such as readings
from various sensors (the observations). The state of the
system evolves from one state (e.g. position and speed at time t)
to another (position and speed at time t+1) according to a known
transition equation, possibly including random perturbations and
intervention effects. The observations are derived from the
state values by a an observation equation (e.g. observation at time t
= position + noise), also possibly including random disturbances and
intervention effects.

The challenge is to obtain the best estimate of the unknown state
considering the set of available observations at a given point in
time. Due to the presence of noise disturbances, it is generally not
possible to simply use the observations directly because they lead to
estimators which are too erratic.
During the 1960s, the Kalman filtering and smoothing algorithm
was developed and popularized to efficiently and optimally solve this
etimation problem.
The technique is based on an iterative procedure in which state values
are successively predicted given the knowledge of the past
observations, and then updated upon the reception of the next
observation. Because of the predict-and-update nature of Kalman
filtering, it can also be interpreted under a Bayesian perspective.
}
\section{Dynamic linear models}{
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

There are two constitutive operations for dynamic linear models:
filtering and smoothing.
In a few words, filtering is the operation consisting in estimating
the state values at time t, using only observations up to (and
including) t-1.
On the contrary, smoothing is the operation which aims at estimating
the state values using the whole set of observations.
}
\details{
The distinguishing aspect of this package is that it provides
functions for naming and extracting components of DLMs (see below),
and a unified interface compatible with other state-space R packages
for filtering and smoothing:
\itemize{
\item package \code{KFAS}: implements exact diffuse initialization and
      supports filtering, smoothing and likelihood computation for the
      exponential family state-space models (as of v0.9.9), used by
      default
\item package \code{dlm}: good general purpose package with many
      helper functions available and some support for Bayesian
      analysis (as of v1.1-2)
\item package \code{FKF}: very fast and memory efficient but has no
      smoothing algorithm (as of v0.1.2)
\item package \code{sspir}: provides (extended) Kalman filter and
      Kalman smoother for models with support for the exponential
      family, but it has no support for the multivariate case,
      exact diffuse initialization or importance sampling,
      and does not support missing values in covariates
      (as of v0.2.8)
}
}
\section{State-space form and notations}{
The state-space model is represented as follows:
\itemize{
\item initial state:        \eqn{alpha(0) \sim N(a(0), P(0))}
\item observation equation: \eqn{y(t)       = Z(t)alpha(t) + eta(t)}
\item transition equation:  \eqn{alpha(t+1) = T(t)alpha(t) + R(t)eps(t)}
\item observation disturbance: \eqn{eta(t) \sim N(0, H(t))}
\item state disturbance:       \eqn{eps(t) \sim N(0, Q(t))}
\item state mean and covariance matrix: \eqn{a(t) = E[alpha(t)]}
and \eqn{P(t) = cov(alpha(t))}
}
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
consisiting of \eqn{level+trend+seasonal+cycle} is an example of
how individual elements can be added together to form a more complete
model. A typical analysis will consider the model as a whole, but also
look into the values of indivdual terms, for example the variations in
the \eqn{level}, and the evolution of the shape of the \eqn{seasonal}
terms. This is also known as seasonal adjustment.

This package introduces a notion called the "component" to facilitate
the analysis of the DLMs. Mathematically speaking, a component is a
named subset of state variables. Components are automatically created
when the model is built, and the package provides functions which
makes it easier to access and analyze their values afterwards:
expectation, variance and prediction/confidence bands.
}
\section{Notes}{
Work is in progress to provide generalized DLM support. 
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
\keyword{ filtering }
\keyword{ smoothing }
\seealso{
Other R packages and functions of interest (in alphabetical order):
\itemize{
\item decompose() from \{stats\}: Decompose a time series into seasonal, trend and irregular components using moving averages. Deals with additive or multiplicative seasonal component.
\item Package \code{dlm}: Maximum likelihood, Kalman filtering and smoothing, and Bayesian analysis of Normal linear State Space models, also known as Dynamic Linear Models.
\item Package \code{dse}: Package dse provides tools for multivariate, linear, time-invariant, time series models. It includes ARMA and state-space representations, and methods for converting between them. It also includes simulation methods and several estimation functions. The package has functions for looking at model roots, stability, and forecasts at different horizons. The ARMA model representaion is general, so that VAR, VARX, ARIMA, ARMAX, ARIMAX can all be considered to be special cases. Kalman filter and smoother estimates can be obtained from the state space model, and state-space model reduction techniques are implemented. An introduction and User's Guide is available in a vignette.
\item Package \code{FKF}: This is a fast and flexible implementation of the Kalman filter, which can deal with NAs. It is entirely written in C and relies fully on linear algebra subroutines contained in BLAS and LAPACK. Due to the speed of the filter, the fitting of high-dimensional linear state space models to large datasets becomes possible. This package also contains a plot function for the visualization of the state vector and graphical diagnostics of the residuals.
\item Package \code{forecast}: Methods and tools for displaying and analysing univariate time series forecasts including exponential smoothing via state space models and automatic ARIMA modelling.
\item HoltWinters() from \{stats\}: Computes Holt-Winters Filtering of a given time series. Unknown parameters are determined by minimizing the squared prediction error.
\item Package \code{KFAS}: Package KFAS provides functions for Kalman filtering, state, disturbance and simulation smoothing, forecasting and simulation of state space models. All functions can use exact diffuse initialisation when distributions of some or all elements of initial state vector are unknown. Filtering, state smoothing and simulation functions use sequential processing algorithm, which is faster than standard approach, and it also allows singularity of prediction error variance matrix. KFAS also contains function for computing the likelihood of exponential family state space models and function for state smoothing of exponential family state space models.
\item Package \code{MARSS}: The MARSS package fits constrained and unconstrained linear multivariate autoregressive state-space (MARSS) models to multivariate time series data.
\item Package \code{sspir}: A glm-like formula language to define dynamic generalized linear models (state space models). Includes functions for Kalman filtering and smoothing. Estimation of variance matrices can be perfomred using the EM algorithm in case of Gaussian models. Read help(sspir) to get started.
\item stl() from \{stats\}: Decompose a time series into seasonal, trend and irregular components using loess, acronym STL.
\item StructTS() and tsSmooth() from \{stats\}: Fit a structural model for a time series by maximum likelihood. Performs fixed-interval smoothing on a univariate time series via a state-space model. Fixed-interval smoothing gives the best estimate of the state at each time point based on the whole observed series.
}

Other software of interest (in alphabetical order):
\itemize{
\item SAS PROC UCM \code{http://www.sas.com/products/ets}
\item SsfPack \code{http://www.ssfpack.com}
\item STAMP \code{http://www.stamp-software.com}
\item S+FinMetrics \code{http://www.insightful.com/products/finmetrics}
\item TRAMO-SEATS \code{http://www.bde.es/servicio/software/econome.htm}
\item X12-ARIMA \code{http://www.census.gov/srd/www/x12a/}
\item X13-ARIMA-SEATS \code{http://www.census.gov/ts/papers/jsm09bcm.pdf}
}
}
\examples{
require(dlmodeler)

# This section illustrates most of the possibilities offered by the
# package by reproducing famous examples of state-space time-series
# analysis found in the litterature.

###############################################
# analysis from Durbin & Koopman book page 32 #
# random walk                                 #
###############################################

# load and show the data
y <- matrix(Nile,nrow=1)
plot(y[1,],type='l')

# y(t)   = a(t) + eta(t)
# a(t+1) = a(t) + eps(t)
# with the parametrization (phi) proposed in the book
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


################################################
# analysis from Durbin & Koopman book page 163 #
# random walk + stochastic seasonal            #
################################################

# load and show the data
y <- matrix(log(Seatbelts[,'drivers']),nrow=1)
plot(y[1,],type='l')

# y(t)    = a(t) + s1(t) + eta(t)
# a(t+1)  = a(t) + eps_L(t)
# s1(t+1) = -s2(t) - s3(t) - ... - s12(t) + eps_S(t)
# s2(t+1) = s1(t)
# s3(t+1) = s2(t), etc.
mod <- dlmodeler.build.structural(
			pol.order=0,
			pol.sigmaQ=NA,
			dseas.order=12,
			dseas.sigmaQ=NA,
			sigmaH=NA,
			name='p163')

# fit the model by maximum likelihood estimation
fit <- dlmodeler.fit(y, mod, method="MLE")

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
# analysis from Durbin & Koopman book page 166          #
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
m1 <- dlmodeler.build.structural(
		pol.order=0,
		dseas.order=12,
		sigmaH=NA, pol.sigmaQ=NA)
m2 <- dlmodeler.build.regression(
		law,
		sigmaQ=0,
		name='law')
m3 <- dlmodeler.build.regression(
		petrolprice,
		sigmaQ=0,
		name='petrolprice')
mod <- dlmodeler.add(dlmodeler.add(m1,m2),m3,name='p166')

# fit the model by maximum likelihood estimation
fit <- dlmodeler.fit(y, mod, method="MLE")

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
mod <- dlmodeler.build.polynomial(0,sigmaQ=NA,name='p32')

# fit the model by maximum likelihood estimation and compute the
# 1-step ahead MSE
fit.mle <- dlmodeler.fit(y, mod, method="MLE")
mean((fit.mle$filtered$f[1,10:100]-y[1,10:100])^2)

# fit the model by minimizing the 1-step ahead MSE
fit.mse1 <- dlmodeler.fit(y, mod, method="MSE", ahead=1, start=10)
mean((fit.mse1$filtered$f[1,10:100]-y[1,10:100])^2)

# fit the model by minimizing the 4-step ahead MSE
fit.mse4 <- dlmodeler.fit(y, mod, method="MSE", ahead=4, start=10)
mean((fit.mse4$filtered$f[1,10:100]-y[1,10:100])^2)

# compare the 1-step ahead forecasts for these models
# as can be expected, the MLE and MSE1 models roughly
# have the same means
plot(y[1,],type='l')
lines(fit.mle$filtered$f[1,],col='dark blue')
lines(fit.mse1$filtered$f[1,],col='dark green')
lines(fit.mse4$filtered$f[1,],col='dark red')


#################################################
# looking at variances and prediction intervals #
#################################################

# load and show the data
y <- matrix(log(Seatbelts[,'drivers']),nrow=1)
plot(y[1,],type='l')

# model with level + seasonal
mod <- dlmodeler.build.structural(
			pol.order=0,
			pol.sigmaQ=NA,
			dseas.order=12,
			dseas.sigmaQ=NA,
			sigmaH=NA,
			name='p163')

# fit the model by maximum likelihood estimation, filter & smooth
fit <- dlmodeler.fit(y, mod, method="MLE")
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
