\name{FStest}
\alias{FStest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
k-Sample FS Test of Equal Distributions
}
\description{
Performs the distribution free exact k-sample test for equality of multivariate distributions in the HDLSS regime.
}
\usage{
FStest(M, labels, sizes, n_clust, randomization = TRUE, clust_alg = "knwClustNo", 
kmax = 2 * n_clust, s_psi = 1, s_h = 1, lb = 1, n_sts = 1000, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
\eqn{n\times d} observations matrix of pooled sample, the observations should be grouped by their respective classes
}
  \item{labels}{
length \eqn{n} vector of membership index of observations
}
  \item{sizes}{
vector of sample sizes
}
  \item{n_clust}{
number of the Populations
}
  \item{randomization}{
logical; if TRUE (default), randomization test and FALSE, non-randomization test
}
  \item{clust_alg}{
\code{"knwClustNo"}(default) or \code{"estclustNo"}(for MFS test); modified K-means algorithm used for clustering
}
  \item{kmax}{
maximum value of total number of clusters to estimate total number of clusters in the whole observations, default: \code{2*n_clust}
}
  \item{s_psi}{
function required for clustering, 1 for \eqn{t^2}, 2 for \eqn{1-\exp(-t)}, 3 for \eqn{1-\exp(-t^2)}, 4 for \eqn{\log(1+t)}, 5 for \eqn{t}
}
  \item{s_h}{
function required for clustering, 1 for \eqn{\sqrt t}, 2 for \eqn{t}
}
  \item{lb}{
each observation is partitioned into some numbers of smaller vectors of same length \eqn{lb}, default: \eqn{1} 
}
  \item{n_sts}{
number of simulation of the test statistic, default: \eqn{1000}
}
  \item{alpha}{
numeric, confidence level \eqn{\alpha}, default: \eqn{0.05}
}
}

\value{
%%  ~Describe the value returned
   FStest returns a list containing the following items:
  \item{estClustLabel}{a vector of length \eqn{n} of estimated class membership index of all observations}
  \item{obsCtyTab}{observed contingency table}
  \item{ObservedProb}{value of the observed test statistic}
  \item{FCutoff}{cut-off of the test}
  \item{randomGamma}{randomized coefficient of the test}
  \item{estPvalue}{estimated p-value of the test}
  \item{decisionF}{if returns \eqn{1}, reject the null hypothesis and if returns \eqn{0}, fails to reject the null hypothesis}
  \item{estClustNo}{total number of the estimated classes}
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
   # muiltivariate normal distribution:
   # generate data with dimension d = 500
   set.seed(151)
   n1=n2=n3=n4=10
   k = 4
   d = 500
   I1 <- matrix(rnorm(n1*d,mean=0,sd=1),n1,d)
   I2 <- matrix(rnorm(n2*d,mean=0.5,sd=1),n2,d) 
   I3 <- matrix(rnorm(n3*d,mean=1,sd=1),n3,d) 
   I4 <- matrix(rnorm(n4*d,mean=1.5,sd=1),n4,d) 
   levels <- c(rep(0,n1), rep(1,n2), rep(2,n3), rep(3,n4)) 
   X <- as.matrix(rbind(I1,I2,I3,I4)) 
   #FS test:
   results <- FStest(M=X, labels=levels, sizes = c(n1,n2,n3,n4), n_clust = k)
  
   ## outputs:
   results$estClustLabel
   #[1] 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3

   results$obsCtyTab
   #      [,1] [,2] [,3] [,4]
   #[1,]   10    0    0    0
   #[2,]    0   10    0    0
   #[3,]    0    0   10    0
   #[4,]    0    0    0   10

   results$ObservedProb
   #[1] 2.125236e-22

   results$FCutoff
   #[1] 1.115958e-07

   results$randomGamma
   #[1] 0

   results$estPvalue
   #[1] 0

   results$decisionF
   #[1] 1

}
\keyword{ FStest function}
