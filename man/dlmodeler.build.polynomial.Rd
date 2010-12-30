\name{dlmodeler.build.polynomial}
\alias{dlmodeler.build.polynomial}
\title{
Build a polynomial model
}
\description{
Builds an univariate polynomial DLM of the specified order.
}
\usage{
dlmodeler.build.polynomial(ord, sigmaH = 1, sigmaQ = 0,
                           name = "polynomial")
}
\arguments{
  \item{ord}{order of the polynomial (\eqn{0} = constant, \eqn{1} = linear, ...).}
  \item{sigmaH}{std dev of the observation disturbance.}
  \item{sigmaQ}{std dev of the state disturbances.}
  \item{name}{an optional name to be given to the resulting DLM.}
}
\details{
The polynomial term is of the form
\eqn{a[1] + a[2]t + a[3]t^2 ... + a[ord]t^ord}.

The initial value \code{P0inf} is parametered to use exact diffuse initialisation
(if supported by the back-end).

Remember: the observation disturbance is generally greater than 0, and the state
disturbance is 0 for classical deterministic models. By looking at the prototype
for this function, you can assume sigmaH refers to the observation and sigmaQ
refers to the state. This rule of thumb is useful to remember how this function
should be used.
}
\value{
An object of class \code{dlmodeler} representing the polynomial model.
}
\references{
Durbin, and Koopman, Time Series Analysis by State Space Methods, Oxford University Press (2001),
pp. 38-45.
}
\author{
Cyrille Szymanski <cnszym@gmail.com>
}
\note{
State representations are generally not unique, so other forms could be used to achieve the same goals.
}

\seealso{
\code{\link{dlmodeler}},
\code{\link{dlmodeler.build}},
\code{\link{dlmodeler.build.dseasonal}},
\code{\link{dlmodeler.build.tseasonal}},
\code{\link{dlmodeler.build.structural}},
\code{\link{dlmodeler.build.arima}},
\code{\link{dlmodeler.build.regression}}
}
\examples{
require(dlmodeler)

# generate some quarterly data
n <- 80
level <- 12
sigma <- .75
season <- c(5,6,8,2)
y <- level + rep(season,n/4) + rnorm(n, mean=0, sd=sigma)

# deterministic level + quarterly seasonal + disturbance
mod1 <- dlmodeler.build.polynomial(0,sigmaH=sigma)
mod2 <- dlmodeler.build.dseasonal(4,sigmaH=0)
mod <- dlmodeler.add(mod1, mod2)
f <- dlmodeler.filter(y, mod)

# show the one step ahead forecasts
plot(y,type='l')
lines(f$f[1,],col='light blue')
}
\keyword{ dlm }
\keyword{ polynomial }