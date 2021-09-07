\name{vamForCollection}
\alias{vamForCollection}
\title{VAM method for multiple gene sets}
\description{
Executes the Variance-adjusted Mahalanobis (VAM) method (\code{\link{vam}}) on multiple gene sets, i.e., a gene set collection.
}
\usage{
    vamForCollection(gene.expr, gene.set.collection, tech.var.prop, 
        gene.weights, center=FALSE, gamma=TRUE)
}
\arguments{
\item{gene.expr}{An n x p matrix of gene expression values for n cells and p genes.}
\item{gene.set.collection}{List of m gene sets for which scores are computed.
    Each element in the list corresponds to a gene set and the list element is a vector
    of indices for the genes in the set. The index value is defined relative to the
    order of genes in the \code{gene.expr} matrix. Gene set names should be specified as list names.}
\item{tech.var.prop}{See description in \code{\link{vam}}} 
\item{gene.weights}{See description in \code{\link{vam}}. 
    If specified as a single vector of weights, weights must be specified for all p genes and
    the same weights are used for all gene sets.
    To use different weights for each set, specify as a list of the same length as the 
    \code{gene.set.collection} list. In this case, each list element should be a vector of 
    gene weights of the same length as the size of the corresponding gene set.}
\item{center}{See description in \code{\link{vam}}} 
\item{gamma}{See description in \code{\link{vam}}} 
}
\value{
A list containing two elements:
    \itemize{
    \item "cdf.value": n x m matrix of 1 minus the one-sided p-values for the m gene sets and n cells. 
   \item "distance.sq": n x m matrix of squared adjusted Mahalanobis distances for the m gene sets and n cells. 
     }
}
\seealso{\code{\link{vam}},\code{\link{vamForSeurat}}} 
\keyword{file}
\examples{
    # Simulate Poisson expression data for 10 genes and 10 cells
    gene.expr=matrix(rpois(100, lambda=2), nrow=10)
    # Simulate technical variance proportions
    tech.var.prop=runif(10)
    # Define a collection with two disjoint sets that span the 10 genes
    collection=list(set1=1:5, set2=6:10)    
    # Execute VAM on both sets using default values for center and gamma
    vamForCollection(gene.expr=gene.expr, gene.set.collection=collection,
        tech.var.prop=tech.var.prop)
    # Create weights that prioritize the first 2 genes for the first set 
    # and the last 2 genes for the second set
    gene.weights = list(c(2,2,1,1,1),c(1,1,1,2,2))
    # Execute VAM using the weights
    vamForCollection(gene.expr=gene.expr, gene.set.collection=collection,
        tech.var.prop=tech.var.prop, gene.weights=gene.weights)
}