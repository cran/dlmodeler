\name{dlmodeler.build}
\alias{dlmodeler.build}
\title{
Build a DLM
}
\description{
Builds a DLM with the supplied design matrices, or an "empty" DLM of the specified dimensions.
}
\usage{
dlmodeler.build(a0 = NULL, P0 = NULL, P0inf = NULL,
                Tt = NULL, Rt = NULL, Qt = NULL,
                Zt = NULL, Ht = NULL,
                dimensions = NULL,
                name = 'noname', components = list())
}
\arguments{
  \item{a0}{initial state vector.}
  \item{P0}{initial state covariance matrix.}
  \item{P0inf}{diffuse part of \code{P0}, matrix of zeros and ones.}
  \item{Tt}{state transition matrix.}
  \item{Rt}{state disturbance selection matrix.}
  \item{Qt}{state disturbance covariance matrix.}
  \item{Zt}{observation design matrix.}
  \item{Ht}{observation disturbance covariance matrix.}
  \item{dimensions}{vector of dimensions \eqn{(m,r,d)}.}
  \item{name}{an optional name to be given to the resulting DLM.}
  \item{components}{optional list of components.}
}
\details{
A DLM can be constructed either by specifying all the elements
\code{a0}, \code{P0}, \code{P0inf},\code{Tt}, \code{Rt}, \code{Qt}, \code{Zt} and \code{Ht}
or by simply giving the \code{dimensions} \eqn{m}, \eqn{r} and \eqn{d}
(in which case the DLM is created with zero-filled elements of the appropriate dimension).

See \code{\link{dlmodeler}} for information about the state-space representation
adopted in this package.

This function is called by the helper functions referenced below.
}
\value{
An object of class \code{dlmodeler} representing the model.
}
\author{
Cyrille Szymanski <cnszym@gmail.com>
}

\seealso{
\code{\link{dlmodeler}}, 
\code{\link{dlmodeler.check}}, 
\code{\link{dlmodeler.build.polynomial}}, 
\code{\link{dlmodeler.build.dseasonal}}, 
\code{\link{dlmodeler.build.tseasonal}}, 
\code{\link{dlmodeler.build.structural}}, 
\code{\link{dlmodeler.build.arima}}, 
\code{\link{dlmodeler.build.regression}}
}
\examples{
\dontrun{
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
}
\keyword{ dlm }
