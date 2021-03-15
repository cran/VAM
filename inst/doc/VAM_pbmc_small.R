### R code from vignette source 'VAM_pbmc_small.Rnw'

###################################################
### code chunk number 1: VAM_pbmc_small.Rnw:15-16
###################################################
library(VAM)


###################################################
### code chunk number 2: VAM_pbmc_small.Rnw:23-31
###################################################
if (requireNamespace("Seurat", quietly=TRUE)) {
	SeuratObject::pbmc_small
	gene.names = rownames(SeuratObject::pbmc_small)
	gene.names[1:5]
	Seurat::VariableFeatures(SeuratObject::pbmc_small)[1:5]
} else {
	message("Seurat package not available! Not executing associated vignette content.")
}	


###################################################
### code chunk number 3: VAM_pbmc_small.Rnw:38-54
###################################################
if (requireNamespace("Seurat", quietly=TRUE)) {
	gene.set.name = "Test"
	gene.ids = c("PPBP", "IGLL5", "VDAC3", "CD1C", "AKR1C3")
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
### code chunk number 4: VAM_pbmc_small.Rnw:61-68
###################################################
if (requireNamespace("Seurat", quietly=TRUE)) {
	pbmc.vam = vamForSeurat(seurat.data=SeuratObject::pbmc_small,
	    gene.set.collection=gene.set.collection,
	    center=F, gamma=T, sample.cov=F, return.dist=T)
} else {
	message("Seurat package not available! Not executing associated vignette content.")
}	


###################################################
### code chunk number 5: VAM_pbmc_small.Rnw:73-79
###################################################
if (requireNamespace("Seurat", quietly=TRUE)) {
	pbmc.vam@assays$VAMdist[1,1:10]
	pbmc.vam@assays$VAMcdf[1,1:10]
} else {
	message("Seurat package not available! Not executing associated vignette content.")
}	


###################################################
### code chunk number 6: VAM_pbmc_small.Rnw:86-97
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


