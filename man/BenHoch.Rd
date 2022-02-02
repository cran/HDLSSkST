\name{BenHoch}
\alias{BenHoch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Benjamini-Hochbergs step-up-procedure (1995)
}
\description{
Benjamini-Hochbergs step-up-procedure (1995) for multiple tests.
}
\usage{
BenHoch(pvalues, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvalues}{
vector of p-values
}
  \item{alpha}{
numeric, false discovery rate controling level \eqn{\alpha}, default: \eqn{0.05}
}
}

\value{
%%  ~Describe the value returned
a vector of \eqn{0}s and \eqn{1}s. \eqn{0}: fails to reject the corresponding hypothesis and \eqn{1}: reject the corresponding hypothesis
}
\references{
Yoav Benjamini and Yosef Hochberg (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing, \emph{Journal of the Royal statistical society: series B (Methodological)} 57.1: 289-300, doi: 10.2307/2346101.
}
\author{ Biplab Paul, Shyamal K. De and Anil K. Ghosh

Maintainer: Biplab Paul<paul.biplab497@gmail.com>
}
\examples{
   # Benjamini-Hochbergs step-up-procedure:
   pvalues <- c(0.50,0.01,0.001,0.69,0.02,0.05,0.0025)
   alpha <- 0.05
   BenHoch(pvalues, alpha)

   ## outputs:
   #[1] 0 1 1 0 1 0 1
}
\keyword{BenHoch function}
