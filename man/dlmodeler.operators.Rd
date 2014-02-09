\name{dlmodeler.operators}
\alias{dlmodeler.add}
\alias{+.dlmodeler}
\alias{dlmodeler.multiply}
\alias{*.dlmodeler}
\alias{dlmodeler.bind}
\alias{\%\%.dlmodeler}
\title{
Add, multiply or bind models
}
\description{
Add two DLMs together, performing an outer sum.

Multiply a model by a numeric constant.

Bind two DLMs together, creating a multi-variate model.
}
\usage{
dlmodeler.add(e1, e2, name = NULL)

\method{+}{dlmodeler}(e1,e2)

dlmodeler.multiply(e1, e2)

\method{*}{dlmodeler}(e1,e2)

dlmodeler.bind(e1, e2, name = NULL)

\method{\%\%}{dlmodeler}(e1, e2)
}
\arguments{
  \item{e1, e2}{an object of class \code{dlmodeler} or a numeric value.}
  \item{name}{an optional name to be given to the resulting DLM.}
}
\details{
Addition: The state vector of the resulting DLM is equal to the concatenation of
the state vectors of \code{mod1} and \code{mod2}. The observation vector of the
resulting DLM is equal to the sum of the observation vectors of \code{mod1} and
\code{mod2}.

Multiplication: the observation vector of the resulting DLM is multiplied
by the supplied numeric constant.

Binding: The state vector of the resulting DLM is equal to the concatenation of
the state vectors of \code{mod1} and \code{mod2}. The observation vector of the
resulting DLM is equal to the concatenation of the observation vectors of
\code{mod1} and \code{mod2}.
}
\value{
An object of class \code{dlmodeler}.
}
\references{
Giovanni Petris, An R Package for Dynamic Linear Models. Journal of Statistical Software, 36(12), 1-16. http://www.jstatsoft.org/v36/i12/.
}
\author{
Cyrille Szymanski <cnszym@gmail.com>
}

\seealso{
\code{\link{dlmodeler}}
}
\examples{
require(dlmodeler)

# create the following model:
# deterministic level + quarterly seasonal + disturbance
mod1 <- dlmodeler.build.polynomial(0,sigmaH=.1) +
        4* dlmodeler.build.dseasonal(4,sigmaH=0)

# create a multivariate model by binding the previous model
# with a stochastic trend model
mod2 <- mod1 \%\% dlmodeler.build.polynomial(1,sigmaH=0)
}
\keyword{ dlm }
