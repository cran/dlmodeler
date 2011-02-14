\name{dlmodeler.add}
\alias{dlmodeler.add}
\title{
Add two models
}
\description{
Add two DLMs together, performing an outer sum.
}
\usage{
dlmodeler.add(mod1, mod2, name = NULL)
}
\arguments{
  \item{mod1, mod2}{objects of class \code{dlmodeler}.}
  \item{name}{an optional name to be given to the resulting DLM.}
}
\details{
The state vector of the resulting DLM is equal to the concatenation of
the state vectors of \code{mod1} and \code{mod2}.

The observation vector of the resulting DLM is equal to the sum of
the observation vectors of \code{mod1} and \code{mod2}.
}
\value{
An object of class \code{dlmodeler} representing the sum of \code{mod1} and \code{mod2}.
}
\references{
Giovanni Petris, An R Package for Dynamic Linear Models. Journal of Statistical Software, 36(12), 1-16. http://www.jstatsoft.org/v36/i12/.
}
\author{
Cyrille Szymanski <cnszym@gmail.com>
}

\seealso{
\code{\link{dlmodeler}}, 
\code{\link{dlmodeler.bind}}
}
\examples{
require(dlmodeler)

# create the following model:
# deterministic level + quarterly seasonal + disturbance
mod1 <- dlmodeler.build.polynomial(0,sigmaH=.1)
mod2 <- dlmodeler.build.dseasonal(4,sigmaH=0)
mod <- dlmodeler.add(mod1, mod2)
}
\keyword{ dlm }
