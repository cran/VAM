#
# VAM.R
#
# Implementation of the Variance Adjusted Mahalanobis (VAM) method which
# computes a modified Mahalanobis distance measure on single cell gene expression data.
#
# @author rob.frost@dartmouth.edu
#

#
# Uses the variance adjusted Mahalanobis (VAM) method to compute distance statistics and one-sided p-values 
# for all cells in the specified gene expression matrix. This matrix should reflect the subset of the full 
# expression profile that corresponds to a single gene set. The p-values will be computed using either a 
# chi-square distribution, a non-central chi-square distribution or gamma distribution as controlled by the
# center and gamma arguments for the one-sided alternative hypothesis that the expression values in the 
# cell are further from the mean (center=T) or origin (center=F) than expected under the null
# of uncorrelated technical noise, i.e., gene expression variance is purely technical and all genes are uncorrelated.
#
# Inputs:
#
# -gene.exprs: An n x p matrix of gene expression values for n cells and p genes. 
# -tech.var.prop: Vector of technical variance proportions for each of the p genes. If specified, the Mahalanobis distance
#          will be computed using a diagonal covariance matrix generated using these proportions. If not specified, the
#          Mahalanobis distances will be computed using a diagonal covariance matrix generated from the sample variances.
# -center: If true will mean center the values in the computation of the Mahalanobis statistic.
#          If false, will compute the Mahalanobis distance from the origin. Default is F.
# -gamma: If true, will fit a gamma distribution to the non-zero squared Mahalanobis distances computed from 
#         a row-permuted version of gene.exprs. The estimated gamma distribution will be used to
#         a one-sided p-value for each sample. If false, will compute the p-value using the standard 
#         chi-square approximation for the squared Mahalanobis distance (or non-central if center=F). Default is T.
#
# Output:
#  
#   A data.frame with the following elements (row names will match row names from gene.expr):
#     -cdf.value: 1 minus the one-sided p-values computed from the squared adjusted Mahalanobis distances.
#     -distance.sq: The squared adjusted Mahalanobis distances for the n cells.
#
#

vam = function(gene.expr, tech.var.prop, center=FALSE, gamma=TRUE){    
  if (missing(gene.expr)) {
    stop("Missing gene expression matrix!")
  }    
    
  #----------------------------------------------------------------------------------------------------  
  # Make sure no genes have a 0 mean value 
  #----------------------------------------------------------------------------------------------------  
  mean.values = apply(gene.expr, 2,  mean)
  zero.mean.genes = which(mean.values == 0)
  if (length(zero.mean.genes) > 0) {
    warning("Removing ", length(zero.mean.genes), " genes with 0 mean values.")
    genes.to.keep = which(mean.values > 0)
    mean.values = mean.values[genes.to.keep]
    gene.expr = gene.expr[,genes.to.keep]
    if (length(genes.to.keep) == 1) {
      # Force vector to matrix
      warning("Gene set has just a single member after removing genes with 0 mean values!")
      gene.expr = as.matrix(gene.expr)
    }            
    if (!missing(tech.var.prop)) {    
      tech.var.prop = tech.var.prop[genes.to.keep]
    }
  }  
    
  #----------------------------------------------------------------------------------------------------  
  # Compute gene variances. Remove any genes with 0 variance
  #----------------------------------------------------------------------------------------------------  
  gene.var = apply(gene.expr, 2, var)
  zero.var.genes = which(gene.var == 0)
  if (length(zero.var.genes) > 0) {
    warning("Removing ", length(zero.var.genes), " genes with 0 variance.")
    genes.to.keep = which(gene.var > 0)
    mean.values = mean.values[genes.to.keep]
    gene.expr = gene.expr[,genes.to.keep]
    if (length(genes.to.keep) == 1) {
      # Force vector to matrix
      warning("Gene set has just a single member after removing genes with 0 variance!")
      gene.expr = as.matrix(gene.expr)
    }        
    gene.var = gene.var[genes.to.keep]
    if (!missing(tech.var.prop)) {        
      tech.var.prop = tech.var.prop[genes.to.keep]
    }
  }
  
  n = nrow(gene.expr)
  p = ncol(gene.expr)    

  #----------------------------------------------------------------------------------------------------    
  # If tech.var.prop was specified, use that to estimate the technical
  # variance, otherwise, just use the sample variance
  #----------------------------------------------------------------------------------------------------    
  if (!missing(tech.var.prop)) {
    if (length(tech.var.prop) != p) {
      stop("Length of tech.var.prop, ", length(tech.var.prop), " not equal to the number of genes, ", p, "!")
    }
    tech.var = tech.var.prop * gene.var 
  } else {
    tech.var = gene.var
  }
  
  #----------------------------------------------------------------------------------------------------
  # Create simple inverse of diagonal covariance matrix using technical variance
  #----------------------------------------------------------------------------------------------------
  if (p == 1) {
    inv.cov = as.matrix(1/tech.var)        
  } else {
    inv.cov = diag(1/tech.var)
  }  
    
#    # Estimate covariance matrix
#    cov.mat = cov(as.matrix(gene.expr))
#    
#    # Since this matrix needs to be inverted, 
#    # check for singularity or near-singularity (is condition number < .Machine$doublet.eps) and adjust 
#    # by adding small delta to diagonal
#    cond.num = rcond(cov.mat)
#    if (is.na(cond.num) | cond.num < .Machine$double.eps) {
#      warn("Computed covariance matrix is near singular, adjusting.")
#      cov.mat =  cov.mat + diag(rep(diag.delta, nrow(cov.mat)))
#    }
#    
#    # Invert the covariance matrix 
#    inv.cov = solve(cov.mat)    
    
  #----------------------------------------------------------------------------------------------------
  # Compute the squared Mahalanobis distance. Use custom logic
  # rather than standard R mahalanobis() to support distances from origin and custom
  # covariance matrix.
  #----------------------------------------------------------------------------------------------------
  mahalanobis.sq = apply(gene.expr, 1, function(r) {
        if (center) {
          return ((r-mean.values) %*% inv.cov %*% (r-mean.values))
        } else {
          return (r %*% inv.cov %*% r)
        }
      })

  # Check if any of the distances are infinte
  if (any(is.infinite(mahalanobis.sq))) {
    stop("Computed ", length(which(is.infinite(mahalanobis.sq))), " infinite distances")
  }  
  
  #----------------------------------------------------------------------------------------------------
  # Compute one-sided p-values using the squared Mahalanobis distances using a
  # chi-square, non-central chi-square or gamma null distribution. 
  #----------------------------------------------------------------------------------------------------
  if (gamma) {

    #------------------------------------------------
    # Compute p-values using a gamma distribution    
    #------------------------------------------------  
        
    # To provide a better null distribution, compute adjusted Mahalanobis distances
    # with sample labels permuted for each gene. This will break any correlation between
    # genes and should remove outlier cells.
    gene.expr.null = apply(gene.expr, 2, function(c) {
          return (sample(c, length(c), replace=F))
        })
    
    # Make sure we remove any covariances (variances and means are not impacted
    # by row permutation)
    #if (missing(tech.var)) {  
    #  inv.cov = solve(diag(diag(cov.mat)))
    #}
    
    # Compute the squared distances on the permuted data 
    mahalanobis.sq.null = apply(gene.expr.null, 1, function(r) {
          if (center) {
            return ((r-mean.values) %*% inv.cov %*% (r-mean.values))
          } else {
            return (r %*% inv.cov %*% r)
          }
        })    

    # Fit a gamma distribution to the non-zero squared distances from the row permuated data
    # using maximum likelihood estimation as implemented by fitdistr() in the MASS package. 
    # Let fitdistr pick initial values (silence NaN warnings).
    nonzero.values = which(mahalanobis.sq.null > 0)
    gamma.fit = try(fitdistr(mahalanobis.sq.null[nonzero.values], "gamma", lower=0.01))
    if (inherits(gamma.fit, "try-error")) {
      warning("Estimation of gamma distribution failed, defaulting p-values to 1")
      p.values = rep(1, length(mahalanobis.sq))
    } else {
      # compute the one-sided p-values
      p.values = pgamma(mahalanobis.sq, shape=gamma.fit$estimate[1], rate=gamma.fit$estimate[2], lower.tail=F)      
    }     

  } else if (center) {
   
    #------------------------------------------------------------------
    # Compute p-values using the standard chi-square null distribution
    #------------------------------------------------------------------
  
    p.values = pchisq(mahalanobis.sq, df=p, lower.tail=F)        
        
  } else {
    
    #------------------------------------------------------------------
    # Compute p-values using a non-central chi-square null distribution
    #------------------------------------------------------------------
  
    # Adjust the mean values used for the ncp for the variance structure
    # to get a non-centrality parameter to use with chi-square distribution
    std.mean.values = sqrt(mean.values %*% inv.cov %*% mean.values)
    ncp = sum(std.mean.values^2)
  
    p.values = pchisq(mahalanobis.sq, df=p, lower.tail=F, ncp=ncp)            
    
  }
    
  # ensure no p-values are 0; set any 0 values to the minimum non-zero value
  zero.pvals = which(p.values == 0)
  if (length(zero.pvals) > 0) {
    min.non.zero = min(p.values[which(p.values > 0)])
    p.values[zero.pvals] = min.non.zero
  }
  
  # Create a data.frame to hold the results
  results = data.frame(cdf.value = 1-p.values,
                       distance.sq = mahalanobis.sq)
  rownames(results) = rownames(gene.expr)                   
  
  return (results)
}

#
# Calls the vam() method for multiple gene sets.
# 
# Inputs:
#
# -gene.expr: A n x p matrix of gene expression values for n cells and p genes.
# -gene.set.collection: List of m gene sets for which scores are computed.
#              Each element in the list corresponds to a gene set and the list element is a vector
#              of indices for the genes in the set. The index value is defined relative to the
#              order of genes in the gene.exprs matrix. Gene set names should be specified as list names.
#              See createGeneSetCollection() for utility function that can be used to 
#              help generate this list of indices.
# -tech.var.prop: See description in vam().
# -center: See description in vam().
# -gamma: See description in vam().
#
# Output:
#  
#   A list containing two elements:
#     -cdf.value: n x m matrix of 1 minus the one-sided p-values for the m gene sets and n samples/cells. 
#     -distance.sq: n x m matrix of squared adjusted Mahalanobis distances for the m gene sets and n samples/cells. 
#
vamForCollection = function(gene.expr, gene.set.collection, tech.var.prop, center=FALSE, gamma=TRUE) {
    
  if (missing(gene.expr)) {
    stop("Missing gene expression matrix!")
  }
  if (missing(gene.set.collection)) {
    stop("Missing gene set collection list!")
  }  
  
  p = ncol(gene.expr)
  n = nrow(gene.expr)
  
  if (!missing(tech.var.prop)) {
    if (length(tech.var.prop) != p) {
      stop("Length of tech.var.prop ", length(tech.var.prop), 
          " does not match the number of genes in the expression matrix ", p)
    }
  }
  
  cell.ids = rownames(gene.expr)
  num.sets = length(gene.set.collection)
  set.names = names(gene.set.collection)
  gene.ids = colnames(gene.expr)
  set.sizes = unlist(lapply(gene.set.collection, length))
  min.set.size = min(set.sizes)
  median.set.size = median(set.sizes)  
    
  # Prepare the result matrices
  results = list()
  results$distance.sq = matrix(0, nrow=n, ncol=num.sets,
      dimnames=list(cell.ids, set.names))
  results$cdf.value = matrix(0, nrow=n, ncol=num.sets,
      dimnames=list(cell.ids, set.names))  
  
  message("Computing VAM distances for ", num.sets, " gene sets, ", n, " cells and ", p, " genes.")
  message("Min set size: ", min.set.size, ", median size: ", median.set.size)

  # Process all gene sets in the collection
  for (i in 1:num.sets) {
    set.members = gene.set.collection[[i]]
    set.size = set.sizes[i]
    if (i %% 50 == 0) {
      message("Computing for gene set ", i, " of size ", set.size)
      #message("Set members: ", paste0(set.members, collapse=","))
    }
    set.exprs = gene.expr[,set.members]   

    if (set.size == 1) {
      # Force vector to matrix
      warning("Gene set ", i, " has just a single member!")
      set.exprs = as.matrix(set.exprs)
    }
    # Execute VAM for this set
    if (!missing(tech.var.prop)) {
      # Work-around that supports a collection object that is a list of vectors of gene IDs vs.
      # list of vectors of gene indices as created via createGeneSetCollection()
      names(tech.var.prop) = colnames(gene.expr)
      vam.results = vam(gene.expr=set.exprs, tech.var.prop=tech.var.prop[set.members],
          center=center, gamma=gamma)
    } else {
      vam.results = vam(gene.expr=set.exprs, center=center, gamma=gamma)
    }
    results$distance.sq[,i] = vam.results$distance.sq
    results$cdf.value[,i] = vam.results$cdf.value
  }    

  return (results)
}

#
# Utility function that creates a gene set collection list in the format required
# by vamForCollection() given the gene IDs measured in the expression matrix and a 
# list of gene sets as defined by the IDs of the member genes.
#
# Inputs:
#
# -gene.ids: Vector of gene IDs. This should correspond to the genes measured in the 
#            gene expression data.
# -gene.set.collection: List of m gene sets where each element in the list corresponds to
#            a gene set and the list element is a vector gene IDs. List names are gene set names.
# -min.size: Minimum gene set size after filtering out genes not in the gene.ids vector. 
#            Gene sets whose post-filtering size is below this are removed from the final
#            collection list. Default is 1 and cannot be set to less than 1.
# -max.size: Maximum gene set size after filtering out genes not in the gene.ids vector. 
#            Gene sets whose post-filtering size is above this are removed from the final
#            collection list. If not specified, no filtering is performed.
#
# Output:
#  
#   Version of the input gene.set.collection list where gene IDs have been replaced by position indices,
#   genes not present in the gene.ids vector have been removed and gene sets failing the min/max size
#   constraints have been removed.
#
createGeneSetCollection = function(gene.ids, gene.set.collection, min.size=1, max.size) {
  
  # min.size must be at least 1
  if (min.size < 1) {
     stop("Invalid min.size value! Must be 1 or greater.")
  }
  # If max size is set, make sure it is not less than min size
  if (!missing(max.size)) {
    if (max.size < min.size) {
      stop("max.size cannot be less than min.size!")
    }          
  }    
  
  num.genes = length(gene.ids)
  if (num.genes < 1) {
    stop("gene.ids must contain at least one genes!")
  }
  
  num.gene.sets = length(gene.set.collection)   
  if (num.gene.sets < 1) {
    stop("gene.set.collection must contain at least one gene set!")
  }  
  
  num.sets = length(gene.set.collection)
  set.names = names(gene.set.collection)
  gene.set.indices = list()
  for (i in 1:num.sets) {
    set.ids = gene.set.collection[[i]]
    # map IDs to indices
    set.indices = unlist(sapply(set.ids, function(x){which(gene.ids == x)}))
    set.size = length(set.indices)
    if (set.size < min.size) {
      next
    }
    if (!missing(max.size)) {
        if (set.size > max.size) {
          next
        }
    }
    current.index = length(gene.set.indices)+1
    gene.set.indices[[current.index]] = set.indices
    names(gene.set.indices)[current.index] = set.names[i]
  }
    
  return (gene.set.indices)
}



