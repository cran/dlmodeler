\name{dlmodeler.bind}
\alias{dlmodeler.bind}
\title{
Bind two models
}
\description{
Bind two DLMs together, creating a multi-variate model.
}
\usage{
dlmodeler.bind(mod1, mod2, name = NULL)
}
\arguments{
  \item{mod1, mod2}{objects of class \code{dlmodeler}.}
  \item{name}{an optional name to be given to the resulting DLM.}
}
\details{
The state vector of the resulting DLM is equal to the concatenation of
the state vectors of \code{mod1} and \code{mod2}.

The observation vector of the resulting DLM is equal to the
concatenation of the observation vectors of \code{mod1} and
\code{mod2}.
}
\value{
An object of class \code{dlmodeler} representing the sum of
\code{mod1} and \code{mod2}.
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
\code{\link{dlmodeler}}, \code{\link{dlmodeler.add}}
}
\examples{
## Example TODO
}
\keyword{ dlm }
