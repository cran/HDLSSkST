\name{RItest}
\alias{RItest}
\title{k-Sample RI Test of Equal Distributions}
\description{
Performs the distribution free exact k-sample test for equality of multivariate distributions in the HDLSS regime.
}
\usage{
RItest(M, labels, sizes, randomization = TRUE, 
clust_alg = "knwClustNo", s_fn = 1, lb = 1, n_sts = 1000, 
alpha = 0.05)
}
\arguments{
\item{M}{\eqn{n\times d} observations matrix of pooled sample, the observations should be grouped by their respective classes}
  \item{labels}{length \eqn{n} vector of membership index of observations}
  \item{sizes}{vector of sample sizes}
  \item{randomization}{logical; if TRUE (default), randomization test and FALSE, non-randomization test}
  \item{clust_alg}{\code{"knwClustNo"}(default) or \code{"estclustNo"}; modified K-means algorithm used for clustering}
  \item{s_fn}{function required for clustering, 1 (default) for \eqn{1-\exp(-t)}, 2 for \eqn{\log(1+t)}, 3 for \eqn{t}}
  \item{lb}{each observation is partitioned into some numbers of smaller vectors of same length \eqn{lb}, default: \eqn{1}  }
  \item{n_sts}{number of simulation of the test statistic, default: \eqn{1000}}
  \item{alpha}{numeric, confidence level \eqn{\alpha}, default: \eqn{0.05}}
}
\value{
%%  ~Describe the value returned
   RItest returns a list containing the following items:
  \item{estClustLabel}{a vector of length \eqn{n} of estimated class membership index of all observations}
  \item{obsCtyTab}{observed contingency table}
  \item{ObservedRI}{value of the observed test statistic}
  \item{RICutoff}{cut-off of the test}
  \item{randomGamma}{randomized coefficient of the test}
  \item{estPvalue}{estimated p-value of the test}
  \item{decisionRI}{if returns \eqn{1}, reject the null hypothesis and if returns \eqn{0}, fails to reject the null hypothesis}
  \item{estClustNo}{total number of the estimated classes}
}
\references{
Biplab Paul, Shyamal K De and Anil K Ghosh (2020).  Distribution-free Exact k-sample Tests for High Dimensional Low Sample Size data Based on Clustering.

William M Rand (1971). Objective criteria for the evaluation of clustering methods, \emph{Journal of the American Statistical association}, 66(336):846-850, doi:10.1080/01621459.1971.10482356.
}
\author{ Biplab Paul, Shyamal K. De and Anil K. Ghosh

Maintainer: Biplab Paul<biplab.paul@niser.ac.in>
}
\examples{
  # muiltivariate normal distribution:
  # generate data with dimension d = 500
  set.seed(151)
  n1=n2=n3=n4=10
  d = 500
  I1 <- matrix(rnorm(n1*d,mean=0,sd=1),n1,d)
  I2 <- matrix(rnorm(n2*d,mean=0.5,sd=1),n2,d) 
  I3 <- matrix(rnorm(n3*d,mean=1,sd=1),n3,d) 
  I4 <- matrix(rnorm(n4*d,mean=1.5,sd=1),n4,d) 
  levels <- c(rep(0,n1), rep(1,n2), rep(2,n3), rep(3,n4)) 
  X <- as.matrix(rbind(I1,I2,I3,I4)) 
  # RI test:
  results <- RItest(M=X, labels=levels, sizes = c(n1,n2,n3,n4))
  
   ## outputs:
   results$estClustLabel
   #[1] 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3

   results$obsCtyTab
   #      [,1] [,2] [,3] [,4]
   #[1,]   10    0    0    0
   #[2,]    0   10    0    0
   #[3,]    0    0   10    0
   #[4,]    0    0    0   10

   results$ObservedRI
   #[1] 0

   results$RICutoff
   #[1] 0.3307692

   results$randomGamma
   #[1] 0

   results$estPvalue
   #[1] 0

   results$decisionRI
   #[1] 1

}
\keyword{ RItest function}