\name{vamForSeurat}
\alias{vamForSeurat}
\title{VAM wrapper for scRNA-seq data processed using the Seurat framework}
\description{
Executes the Variance-adjusted Mahalanobis (VAM) method (\code{\link{vamForCollection}}) on 
    normalized scRNA-seq data stored in a Seurat object.
    If the Seurat \code{NormalizeData} method was used for normalization, the technical variance of each gene is computed as
    the proportion of technical variance (from \code{FindVariableFeatures}) multiplied by the variance of the  normalized counts. 
    If \code{SCTransform} was used for normalization, the technical variance for each gene is set
    to 1 (the normalized counts output by \code{SCTransform} should have variance 1 if there is only technical variation).
}
\usage{
    vamForSeurat(seurat.data, gene.set.collection, center=FALSE, gamma=TRUE, 
        sample.cov=FALSE, return.dist=FALSE)
}
\arguments{
\item{seurat.data}{The Seurat object that holds the scRNA-seq data. Assumes normalization has already been performed.}
\item{gene.set.collection}{List of m gene sets for which scores are computed.
    Each element in the list corresponds to a gene set and the list element is a vector
    of indices for the genes in the set. The index value is defined relative to the
    order of genes in the relevant \code{seurat.data} Assay object. Gene set names should be specified as list names.}
\item{center}{See description in \code{\link{vam}}} 
\item{gamma}{See description in \code{\link{vam}}} 
\item{sample.cov}{If true, will use the a diagonal covariance matrix generated from the 
    sample variances to compute the squared adjusted Mahalanobis distances (this is equivalent to not specifying
    \code{tech.var} for the \code{\link{vam}} method). If false (default), will use the technical variances
    as determined based on the type of Seurat normalization.}
\item{return.dist}{If true, will return the squared adjusted Mahalanobis distances in a new Assay object called "VAM.dist". 
    Default is F.}
}
\value{
Updated Seurat object that hold the VAM results in one or two new Assay objects:        
    \itemize{
    \item If \code{return.dist} is true, the matrix of squared adjusted Mahalanobis distances will be stored in new 
    Assay object called "VAM.dist". 
    \item The matrix of CDF values (1 minus the one-sided p-values) will be stored in new Assay object called "VAM.cdf". 
     }
}
\seealso{\code{\link{vam}},\code{\link{vamForCollection}}} 
\keyword{file}
\examples{
    # Only run example code if Seurat package is available
    if (requireNamespace("Seurat", quietly=TRUE)) {
        # Define a collection with one gene set for the first 10 genes
        collection=list(set1=1:10)
        # Execute on the pbmc_small scRNA-seq data set included with Seurat
        # See vignettes for more detailed Seurat examples
        vamForSeurat(seurat.data=Seurat::pbmc_small,
            gene.set.collection=collection)
    }
}