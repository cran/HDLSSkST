\name{Holm}
\alias{Holm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Holm's step-down-procedure (1979)
}
\description{
Holm's step-down-procedure (1979) for mutiple tests.
}
\usage{
Holm(pvalues, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvalues}{
vector of p-values
}
  \item{alpha}{
numeric, family wise error rate controling level \eqn{\alpha}, default: \eqn{0.05}
}
}

\value{
%%  ~Describe the value returned
a vector of \eqn{0}s and \eqn{1}s. \eqn{0}: fails to reject the corresponding hypothesis and \eqn{1}: reject the corresponding hypothesis
}
\references{
Sture Holm (1979). A simple sequentially rejective multiple test procedure, \emph{Scandinavian journal of statistics}, 65-70, doi:10.2307/4615733.
}
\author{ Biplab Paul, Shyamal K. De and Anil K. Ghosh

Maintainer: Biplab Paul<biplab.paul@niser.ac.in>
}
\examples{
   # Holm's step down procedure:
   pvalues <- c(0.50,0.01,0.001,0.69,0.02,0.05,0.0025)
   alpha <- 0.05
   Holm(pvalues, alpha)

   ## outputs:
   #[1] 0 0 1 0 0 0 1
   
}
\keyword{Holm function}
