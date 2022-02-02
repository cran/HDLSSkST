\name{AFStest}
\alias{AFStest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
k-Sample AFS Test of Equal Distributions
}
\description{
Performs the distribution free exact k-sample test for equality of multivariate distributions in the HDLSS regime. This an aggregate test of the two sample versions of the FS test over \eqn{\frac{k(k-1)}{2}} numbers of two-sample comparisons, and the test statistic is the minimum of these two sample FS test statistics. Holm's step-down-procedure (1979) and Benjamini-Hochberg procedure (1995) are applied for multiple testing.
}
\usage{
AFStest(M, sizes, randomization = TRUE, clust_alg = "knwClustNo", kmax = 4,
multTest = "Holm", s_psi = 1, s_h = 1, lb = 1, n_sts = 1000, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
\eqn{n\times d} observations matrix of pooled sample, the observations should be grouped by their respective classes
}
  \item{sizes}{
vector of sample sizes
}
  \item{randomization}{
logical; if TRUE (default), randomization test and FALSE, non-randomization test
}
  \item{clust_alg}{
\code{"knwClustNo"}(default) or \code{"estclustNo"}; modified K-means algorithm used for clustering
}
  \item{kmax}{
maximum value of total number of clusters to estimate total number of clusters for two-sample comparition, default: \code{4}
}
  \item{multTest}{
\code{"HOlm"}(default) or \code{"BenHoch"}; different multiple tests
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
   AFStest returns a list containing the following items:
  \item{AFSStat}{value of the observed test statistic}
  \item{AFCutoff}{cut-off of the test}
  \item{randomGamma}{randomized coefficient of the test}
  \item{decisionAFS}{if returns \eqn{1}, reject the null hypothesis and if returns \eqn{0}, fails to reject the null hypothesis}
  \item{multipleTest}{indicates where two populations are different according to multiple tests}

}
\references{
Biplab Paul, Shyamal K De and Anil K Ghosh (2021).  Some clustering based exact distribution-free k-sample tests
applicable to high dimension, low sample size data, \emph{Journal of Multivariate Analysis}, doi:10.1016/j.jmva.2021.104897.

Cyrus R Mehta and Nitin R Patel (1983). A network algorithm for performing Fisher's exact test in rxc contingency tables, \emph{Journal of the American Statistical Association}, 78(382):427-434, doi:10.2307/2288652.

Sture Holm (1979). A simple sequentially rejective multiple test procedure, \emph{Scandinavian journal of statistics}, 65-70, doi:10.2307/4615733.

Yoav Benjamini and Yosef Hochberg (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing, \emph{Journal of the Royal statistical society: series B (Methodological)} 57.1: 289-300, doi: 10.2307/2346101.
}
\author{ Biplab Paul, Shyamal K. De and Anil K. Ghosh

Maintainer: Biplab Paul<paul.biplab497@gmail.com>
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
  X <- as.matrix(rbind(I1,I2,I3,I4)) 
  #AFS test:
  results <- AFStest(M=X, sizes = c(n1,n2,n3,n4))
  
   ## outputs:
   results$AFSStat
   #[1] 5.412544e-06

   results$AFCutoff
   #[1] 0.0109604

   results$randomGamma
   #[1] 0

   results$decisionAFS
   #[1] 1

   results$multipleTest
   #  Population.1 Population.2 rejected pvalues
   #1            1            2     TRUE       0
   #2            1            3     TRUE       0
   #3            1            4     TRUE       0
   #4            2            3     TRUE       0
   #5            2            4     TRUE       0
   #6            3            4     TRUE       0

}
\keyword{ AFStest function}