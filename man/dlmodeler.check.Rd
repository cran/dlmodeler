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
\item{1}{an integer indicating whether the model is valid (1) or not (0)}
\item{2}{dimension of state vector \eqn{m}}
\item{2}{dimension of state disturbance covariance matrix \eqn{r}}
\item{3}{dimension of observation vector \eqn{d}}
\item{4}{an integer indicating if the model has time-varying terms (1) or not (0)}
\item{5}{the number of time steps in Tt, or NA if the matrix is constant}
\item{6}{the number of time steps in Rt, or NA if the matrix is constant}
\item{7}{the number of time steps in Qt, or NA if the matrix is constant}
\item{8}{the number of time steps in Zt, or NA if the matrix is constant}
\item{8}{the number of time steps in Ht, or NA if the matrix is constant}
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
\code{\link{dlmodeler}}, \code{\link{dlmodeler.build}}
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
dlmodeler.check(mod)[1]==1

# an empty DLM with 4 state variables (3 of which are stocastic)
# and bi-variate observations
mod <- dlmodeler.build(dimensions=c(4,3,2))
# print the model
mod
# check if it is valid
dlmodeler.check(mod)[1]==1
}
\keyword{ dlm }
