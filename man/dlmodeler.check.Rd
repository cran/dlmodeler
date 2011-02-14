\name{dlmodeler.check}
\alias{dlmodeler.check}
\title{
Check dimensions and validity
}
\description{
Checks a \code{dlmodeler} object, in particular
for the consistency of the dimensions of its elements.
}
\usage{
dlmodeler.check(model, yt = NULL)
}
\arguments{
  \item{model}{an object of class \code{dlmodeler} to check.}
  \item{yt}{an optional data vector to check with the model.}
}
\details{
See \code{\link{dlmodeler}} for information about the state-space representation
adopted in this package.
}
\value{
An array with the following information:
\item{status}{a boolean indicating whether the model is valid or not}
\item{m}{dimension of state vector \eqn{m}}
\item{r}{dimension of state disturbance covariance matrix \eqn{r}}
\item{d}{dimension of observation vector \eqn{d}}
\item{timevar}{a boolean indicating if the model has time-varying terms or not}
\item{timevar.Tt}{the number of time steps in Tt, or NA if the matrix is constant}
\item{timevar.Rt}{the number of time steps in Rt, or NA if the matrix is constant}
\item{timevar.Qt}{the number of time steps in Qt, or NA if the matrix is constant}
\item{timevar.Zt}{the number of time steps in Zt, or NA if the matrix is constant}
\item{timevar.Ht}{the number of time steps in Ht, or NA if the matrix is constant}
}
\author{
Cyrille Szymanski <cnszym@gmail.com>
}

\seealso{
\code{\link{dlmodeler}}, 
\code{\link{dlmodeler.build}}
}
\examples{
require(dlmodeler)

# a stochastic level+trend DLM
mod <- dlmodeler.build(
		a0 = c(0,0), # initial state: (level, trend)
		P0 = diag(c(0,0)), # initial state variance set to...
		P0inf = diag(2), # ...use exact diffuse initialization
		matrix(c(1,0,1,1),2,2), # state transition matrix
		diag(c(1,1)), # state disturbance selection matrix
		diag(c(.5,.05)), # state disturbance variance matrix
		matrix(c(1,0),1,2), # observation design matrix
		matrix(1,1,1) # observation disturbance variance matrix
)
# print the model
mod
# check if it is valid
dlmodeler.check(mod)$status

# an empty DLM with 4 state variables (3 of which are stocastic)
# and bi-variate observations
mod <- dlmodeler.build(dimensions=c(4,3,2))
# print the model
mod
# check if it is valid
dlmodeler.check(mod)$status
}
\keyword{ dlm }
