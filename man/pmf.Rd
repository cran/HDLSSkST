\name{pmf}
\alias{pmf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized Hypergeometric Probability
}
\description{
A function that provides the probability of observing an \eqn{r\times c} contingency table using generalized hypergeometric probability.
}
\usage{
pmf(M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
\eqn{r\times c} contingency table
}
}

\value{
%%  ~Describe the value returned
a single value between \eqn{0} and \eqn{1}
}
\references{
Biplab Paul, Shyamal K De and Anil K Ghosh (2021).  Some clustering based exact distribution-free k-sample tests
applicable to high dimension, low sample size data, \emph{Journal of Multivariate Analysis}, doi:10.1016/j.jmva.2021.104897.

Cyrus R Mehta and Nitin R Patel (1983). A network algorithm for performing Fisher's exact test in rxc contingency tables, \emph{Journal of the American Statistical Association}, 78(382):427-434, doi:10.2307/2288652.
}
\author{ Biplab Paul, Shyamal K. De and Anil K. Ghosh

Maintainer: Biplab Paul<paul.biplab497@gmail.com>
}
\examples{
   # Generalized hypergeometric probability of rxc Contingency Table:
   mat <- matrix(1:20,5,4, byrow = TRUE)
   pmf(mat)

   ## outputs:
   #[1] 4.556478e-09
}
\keyword{ pmf function}