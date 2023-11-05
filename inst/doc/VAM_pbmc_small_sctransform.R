### R code from vignette source 'VAM_pbmc_small_sctransform.Rnw'

###################################################
### code chunk number 1: VAM_pbmc_small_sctransform.Rnw:15-16
###################################################
library(VAM)


###################################################
### code chunk number 2: VAM_pbmc_small_sctransform.Rnw:23-30
###################################################
if (requireNamespace("Seurat", quietly=TRUE)) {
	SeuratObject::pbmc_small
	gene.names = rownames(SeuratObject::pbmc_small)
	gene.names[1:5]
} else {
	message("Seurat package not available! Not executing associated vignette content.")
}


###################################################
### code chunk number 3: VAM_pbmc_small_sctransform.Rnw:35-44
###################################################
if (requireNamespace("Seurat", quietly=TRUE)) {
	pbmc_sctransform = Seurat::SCTransform(SeuratObject::pbmc_small, verbose=F)
	# Compute PCA and UMAP on the normalized values
	pbmc_sctransform = Seurat::RunPCA(pbmc_sctransform, npcs=10)
	pbmc_sctransform = Seurat::RunUMAP(pbmc_sctransform, dims = 1:10) 
	Seurat::VariableFeatures(pbmc_sctransform)[1:5]
} else {
	message("Seurat package not available! Not executing associated vignette content.")
}


###################################################
### code chunk number 4: VAM_pbmc_small_sctransform.Rnw:51-67
###################################################
if (requireNamespace("Seurat", quietly=TRUE)) {
	gene.set.name = "Test"
	gene.ids = c("NKG7", "PPBP", "GNLY", "PF4", "GNG11")
	# Create a collection list for this gene set
	gene.set.id.list = list()
	gene.set.id.list[[1]] = gene.ids
	names(gene.set.id.list)[1] = gene.set.name
	gene.set.id.list
	# Create the list of gene indices required by vamForSeurat()
	(gene.set.collection = createGeneSetCollection(gene.ids=gene.names,
		gene.set.collection=gene.set.id.list))
	gene.indices = gene.set.collection[[1]]
	(gene.names = gene.names[gene.indices])
} else {
	message("Seurat package not available! Not executing associated vignette content.")
}


###################################################
### code chunk number 5: VAM_pbmc_small_sctransform.Rnw:74-81
###################################################
if (requireNamespace("Seurat", quietly=TRUE)) {
	pbmc.vam = vamForSeurat(seurat.data=pbmc_sctransform,
	    gene.set.collection=gene.set.collection,
	    center=F, gamma=T, sample.cov=F, return.dist=T)
} else {
	message("Seurat package not available! Not executing associated vignette content.")
}


###################################################
### code chunk number 6: VAM_pbmc_small_sctransform.Rnw:86-92
###################################################
if (requireNamespace("Seurat", quietly=TRUE)) {
	pbmc.vam@assays$VAMdist@data[1,1:10]
	pbmc.vam@assays$VAMcdf@data[1,1:10]
} else {
	message("Seurat package not available! Not executing associated vignette content.")
}


###################################################
### code chunk number 7: VAM_pbmc_small_sctransform.Rnw:99-110
###################################################
if (requireNamespace("Seurat", quietly=TRUE)) {
	Seurat::DefaultAssay(object = pbmc.vam) = "VAMcdf"
	Seurat::FeaturePlot(pbmc.vam, reduction="tsne", features=gene.set.name)
} else {
	message("Seurat package not available! Not executing associated vignette content.")
		par(mar = c(0,0,0,0))
	plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
	text(x = 0.5, y = 0.5,paste("Seurat package not available!\n",
					 "FeaturePlot not generated."),
	cex = 1.6, col = "black")
}


