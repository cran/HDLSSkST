\name{HDLSSkST-package}
\alias{HDLSSkST-package}
\docType{package}
\title{
    Distribution-Free Exact High Dimensional Low Sample Size k-Sample Tests
}
\description{
Testing homogeneity of k (\eqn{\geq 2}) multivariate distributions is a classical and challenging problem in statistics, and this becomes even more challenging when the dimension of the data exceeds the sample size. We construct some tests for this purpose which are exact level (size) \eqn{\alpha} tests based on clustering. These tests are easy to implement and distribution-free in finite sample situations. Under appropriate regularity conditions, these tests have the consistency property in HDLSS asymptotic regime, where the dimension of data \eqn{d} grows to \eqn{\infty} while the sample size remains fixed. We also consider a multiscale approach, where the results for the different number of partitions are aggregated judiciously. This package includes eight tests, namely (i) RI test, (ii) FS test, (iii) MRI test, (iv) MFS test, (v) MTRI test , (vi) MTFS test, (vii) ARI test and (viii) AFS test. In MRI and MFS test, we modified the RI and FS test, respectively, using an estimated clustering number. In the multiscale approach (MTRI and MTFS), we use Holm's step-down-procedure (1979) and Benjamini-Hochberg FDR controlling procedure (1995).
}
\author{
Biplab Paul, Shyamal K. De and Anil K. Ghosh

Maintainer: Biplab Paul<biplab.paul@niser.ac.in>
}
\references{
Biplab Paul, Shyamal K De and Anil K Ghosh (2020).  Some clustering based exact distribution-free k-sample tests
applicable to high dimension, low sample size data.

Soham Sarkar and Anil K Ghosh (2019). On perfect clustering of high dimension, low sample size data, \emph{IEEE transactions on pattern analysis and machine intelligence}, doi:10.1109/TPAMI.2019.2912599.

William M Rand (1971). Objective criteria for the evaluation of clustering methods, \emph{Journal of the American Statistical association}, 66(336):846-850, doi:10.1080/01621459.1971.10482356.

Cyrus R Mehta and Nitin R Patel (1983). A network algorithm for performing Fisher's exact test in rxc contingency tables, \emph{Journal of the American Statistical Association}, 78(382):427-434, doi:10.2307/2288652.

Joseph C Dunn (1973). A fuzzy relative of the isodata process and its use in detecting compact well-separated clusters, doi:10.1080/01969727308546046.

Sture Holm (1979). A simple sequentially rejective multiple test procedure, \emph{Scandinavian journal of statistics}, 65-70, doi:10.2307/4615733.

Yoav Benjamini and Yosef Hochberg (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing, \emph{Journal of the Royal statistical society: series B (Methodological)} 57.1: 289-300, doi: 10.2307/2346101.
}
\keyword{ package }
\keyword{ HDLSS}
