\name{print.dlmodeler}
\alias{print.dlmodeler}
\title{
Print a model
}
\description{
Prints a short message describing a DLM.
}
\usage{
\method{print}{dlmodeler}(x,...)
}
\arguments{
  \item{x}{model to be printed.}
  \item{...}{unused.}
}
\details{
This function will print the dimensions of the DLM, the time-varying terms, and the
names of the components.
}
\value{
No value.
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
\keyword{ print }
