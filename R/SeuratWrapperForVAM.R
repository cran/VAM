#
# SeuratWrapperForVAM.R
#
# A Seurat-specific wrapper around VAM.R.
#
# @author rob.frost@dartmouth.edu
#

#
# Computes the adjusted Mahalanobis distance measures using the normalized gene expression values for the cells in a Seurat object.
# If NormalizeData() was used for normalization, the technical variance of each gene is computed as
# the proportion of technical variance (from FindVariableFeatures()) multiplied by the variance of the 
# normalized counts. If SCTransform was used for normalization, the technical variance for each gene is set
# to 1 (the normalized counts output by SCTransform should have variance 1 if there is only technical variation).
# 
# Inputs:
#
# -seurat.data: The Seurat object that holds the scRNA-seq data. 
#  Assumes normalization has already been performed.
# -gene.set.collection: List of m gene sets for which scores are computed.
#  Each element in the list corresponds to a gene set and the list element is a vector
#  of indices for the genes in the set. The index value is defined relative to the
#  order of genes in the gene.exprs matrix. Gene set names should be specified as list names.
#  See createGeneSetCollection() function in VAM.R for utility function that can be used to 
#  help generate this list of indices.
# -center: See vam() function in VAM.R for details.
# -gamma: See vam() function in VAM.R for details. 
# -sample.cov: If true, will use the sample covariance matrix in the Mahalanobis computation. If false (default),
#            will use just the technical variance.
# -return.dist: If true, will return the squared adjusted Mahalanobis distances in a new Assay object called "VAM.dist".
#
# Output:
#
#   Updated Seurat object. 
#   If return.dist is true, the matrix of squared adjusted Mahalanobis distances will be stored in new Assay object called "VAM.dist".
#   The matrix of CDF values (1 minus the one-sided p-values) will be stored in new Assay object called "VAM.cdf".
#
#
vamForSeurat = function(seurat.data,
    gene.set.collection, 
    center=FALSE,
    gamma=TRUE,
    sample.cov=FALSE,
    return.dist=FALSE
) {
  
  if (!requireNamespace("Seurat", quietly=TRUE)) {
    stop("Seurat package not available!")
  }
  
  if (missing(seurat.data)) {
    stop("Missing Seurat object!")
  }
  
  if (missing(gene.set.collection)) {
    stop("Missing gene set collection list!")
  }
  
  # Use the active.assay slot to determine whether SCTransform or the standard normalization logic was used
  if (seurat.data@active.assay == "RNA") {

    # Use the log normalized counts
    normalized.counts = seurat.data@assays$RNA@data
    
    if (!sample.cov) {
      # For the standard normalization pipeline, 
      # get the proportion of variance that appears technical according to FindVariableFeatures
      tech.var.prop = getTechVarProp(seurat.data)
    }
    
  } else if (seurat.data@active.assay == "SCT") {
    
    # Use the log1p transformed corrected counts. The SCT correction process reverses the regression model
    # to generate counts that approximate what would be found if all cells were sequenced to the same depth.
    normalized.counts = seurat.data@assays$SCT@data    
    
    if (!sample.cov) {
      # Get the proportion of variance that appears technical according to the Pearson residuals (these
      # will have variance 1 if entirely technical)
      tech.var.prop = getTechVarPropForSCT(seurat.data)
    }
        
  } else {    
    stop("Unsupported active assay: ", seurat.data@active.assay)
  }    
    
  # Compute the gene set-specific and variance-adjusted Mahalanobis matrix on the normalized counts
  if (sample.cov) {
    vam.results = vamForCollection(gene.expr=Matrix::t(normalized.counts),
        gene.set.collection=gene.set.collection,
        center=center, gamma=gamma)  
  } else {
    vam.results = vamForCollection(gene.expr=Matrix::t(normalized.counts),
        gene.set.collection=gene.set.collection,
        tech.var.prop = tech.var.prop, center=center, gamma=gamma)
  }
  
  # Create Assay object to store the squared distances
  if (return.dist) {
    dist.assay = Seurat::CreateAssayObject(counts = t(vam.results$distance.sq))
    seurat.data[["VAMdist"]] = dist.assay
  }
  
  # Create Assay object to store the CDF values
  cdf.assay = Seurat::CreateAssayObject(counts = t(vam.results$cdf.value))
  seurat.data[["VAMcdf"]] = cdf.assay    
  
  return (seurat.data)
}  


#
# Helper method that computes proportion of technical variance for each gene
# based on the standard normalization method.
#
getTechVarProp = function(seurat.data) {
  
  p = nrow(seurat.data@assays$RNA@data)
  
  # Ensure we have vst results
  meta.features = colnames(seurat.data@assays$RNA@meta.features)
  if (length(meta.features) <= 3 
      | meta.features[2] != "vst.variance" 
      | meta.features[3] != "vst.variance.expected") { 
    message("Did not find vst variance decomposition, setting technical variance proportion to 1")
    return (rep(1, p))
  }
  
  # For the standard normalization pipeline, 
  # get the proportion of variance that appears technical according to FindVariableFeatures
  tech.var.prop = seurat.data@assays$RNA@meta.features[,3]/seurat.data@assays$RNA@meta.features[,2]		
  
  # For NaN entries (vst var is 0), set prop.tech.var to 1
  tech.var.prop[which(is.nan(tech.var.prop))] = 1  
    
  return (tech.var.prop)  
}

#
# Helper method that computes the proportion of technical variance for each gene
# using the SCTransform Pearson residuals
#
getTechVarPropForSCT = function(seurat.data) {
  # make sure that vst.out exists (it 
  sct.misc.elements = names(seurat.data@assays$SCT@misc)
  if (!("vst.out" %in% sct.misc.elements)) {
    stop("Missing seurat.data@assays$SCT@misc$vst.out!",
        "This may be due to the use of the Seurat integration framework (not currently ",
        "supported by VAM) to combine multiple scRNA-seq data sets. For VAM analysis of ",
        "multiple scRNA-seq data sets, SCTransform should be called once on the Seurat object ",
        "generated by merge().")
  }
  
  # The Pearson residuals should have 1 variance if just noise so can use the 
  # saved residual variance to compute the proportion of technical variance
  tech.var.prop = 1/seurat.data@assays$SCT@misc$vst.out$gene_attr$residual_variance
  
  return (tech.var.prop)  
}



