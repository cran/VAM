\name{vam}
\alias{vam}
\title{Variance-adjusted Mahalanobis (VAM) algorithm}
\description{
Implementation of the Variance-adjusted Mahalanobis (VAM) method, which computes distance statistics and one-sided p-values 
for all cells in the specified single cell gene expression matrix. This matrix should reflect the subset of the full 
expression profile that corresponds to a single gene set. The p-values will be computed using either a 
chi-square distribution, a non-central chi-square distribution or gamma distribution as controlled by the
\code{center} and \code{gamma} arguments for the one-sided alternative hypothesis that the expression values in the 
cell are further from the mean (\code{center=T}) or origin (\code{center=F}) than expected under the null 
of uncorrelated technical noise, i.e., gene expression variance is purely technical and all genes are uncorrelated.
}
\usage{
    vam(gene.expr, tech.var.prop, center=FALSE, gamma=TRUE)
}
\arguments{
\item{gene.expr}{An n x p matrix of gene expression values for n cells and p genes.} 
\item{tech.var.prop}{Vector of technical variance proportions for each of the p genes. 
    If specified, the Mahalanobis distance will be computed using a diagonal covariance matrix generated using these proportions. 
    If not specified, the Mahalanobis distances will be computed using a diagonal covariance matrix generated from the 
    sample variances.}
\item{center}{If true will mean center the values in the computation of the Mahalanobis statistic.
    If false, will compute the Mahalanobis distance from the origin. Default is F.}
\item{gamma}{If true, will fit a gamma distribution to the non-zero squared Mahalanobis distances computed from 
    a row-permuted version of \code{gene.expr}. The estimated gamma distribution will be used to compute
    a one-sided p-value for each cell. If false, will compute the p-value using the standard 
    chi-square approximation for the squared Mahalanobis distance (or non-central if \code{center=F}). Default is T.}
}
\value{
A \code{data.frame} with the following elements (row names will match row names from gene.expr):
    \itemize{
      \item "cdf.value": 1 minus the one-sided p-values computed from the squared adjusted Mahalanobis distances. 
      \item "distance.sq": The squared adjusted Mahalanobis distances for the n cells. 
     }
}
\seealso{\code{\link{vamForCollection}},\code{\link{vamForSeurat}}} 
\keyword{file}
\examples{
    # Simulate Poisson expression data for 10 genes and 10 cells
    gene.expr=matrix(rpois(100, lambda=2), nrow=10)
    # Simulate technical variance proportions
    tech.var.prop=runif(10)
    # Execute VAM to compute scores for the 10 genes on each cell
    vam(gene.expr=gene.expr, tech.var.prop=tech.var.prop)
}
