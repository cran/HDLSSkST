\name{MTRItest}
\alias{MTRItest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
k-Sample MTRI Test of Equal Distributions
}
\description{
Performs the distribution free exact k-sample test for equality of multivariate distributions in the HDLSS regime. This test is a multiscale approach based on RI test, where the results for different number of partitions are aggregated judiciously.
}
\usage{
MTRItest(M, labels, sizes, k_max, multTest = "Holm", s_psi = 1, s_h = 1, 
lb = 1, n_sts = 1000, alpha = 0.05)
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
  \item{k_max}{
maximum value of total number of clusters which is required for the test
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
   MTRItest returns a list containing the following items:
  \item{RIvec}{a vector of the Rand indices based on different number of clusters}
  \item{Pvalues}{a vector of RI test p-values based on different number of clusters}
  \item{decisionMTRI }{if returns \eqn{1}, reject the null hypothesis and if returns \eqn{0}, fails to reject the null hypothesis}
  \item{contTabs}{a list of the observed contingency table based on different number of clusters}
  \item{mulTestdec}{a vector of \eqn{0}s and \eqn{1}s. \eqn{0}: fails to reject the corresponding hypothesis and \eqn{1}: reject the corresponding hypothesis}
}
\references{
%% ~put references to the literature/web site here ~
Biplab Paul, Shyamal K De and Anil K Ghosh (2021).  Some clustering based exact distribution-free k-sample tests
applicable to high dimension, low sample size data, \emph{Journal of Multivariate Analysis}, doi:10.1016/j.jmva.2021.104897.

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
  levels <- c(rep(0,n1), rep(1,n2), rep(2,n3), rep(3,n4)) 
  X <- as.matrix(rbind(I1,I2,I3,I4)) 
  #MTRI test:
  results <- MTRItest(X, levels, c(n1,n2,n3,n4), 8)
  
  ## outputs:
  results$RIvec
  #[1] 0.25641026 0.14871795 0.00000000 0.03076923 0.05128205 0.08333333 0.10384615

  results$Pvalues
  #[1] 0 0 0 0 0 0 0

  results$decisionMTRI
  #[1] 1

  results$contTabs
  #$contTabs[[1]]
  #     [,1] [,2]
  #[1,]   10    0
  #[2,]   10    0
  #[3,]    0   10
  #[4,]    0   10

  #$contTabs[[2]]
  #     [,1] [,2] [,3]
  #[1,]   10    0    0
  #[2,]    0   10    0
  #[3,]    0   10    0
  #[4,]    0    0   10

  #$contTabs[[3]]
  #     [,1] [,2] [,3] [,4]
  #[1,]   10    0    0    0
  #[2,]    0   10    0    0
  #[3,]    0    0   10    0
  #[4,]    0    0    0   10

  #$contTabs[[4]]
  #     [,1] [,2] [,3] [,4] [,5]
  #[1,]   10    0    0    0    0
  #[2,]    0   10    0    0    0
  #[3,]    0    0    4    6    0
  #[4,]    0    0    0    0   10

  #$contTabs[[5]]
  #     [,1] [,2] [,3] [,4] [,5] [,6]
  #[1,]   10    0    0    0    0    0
  #[2,]    0   10    0    0    0    0
  #[3,]    0    0    4    6    0    0
  #[4,]    0    0    0    0    8    2

  #$contTabs[[6]]
  #     [,1] [,2] [,3] [,4] [,5] [,6] [,7]
  #[1,]   10    0    0    0    0    0    0
  #[2,]    0    5    5    0    0    0    0
  #[3,]    0    0    0    4    6    0    0
  #[4,]    0    0    0    0    0    8    2

  #$contTabs[[7]]
  #     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
  #[1,]    8    2    0    0    0    0    0    0
  #[2,]    0    0    5    5    0    0    0    0
  #[3,]    0    0    0    0    4    6    0    0
  #[4,]    0    0    0    0    0    0    8    2


  results$mulTestdec
  #[1] 1 1 1 1 1 1 1
}
\keyword{ MTRItest function}
