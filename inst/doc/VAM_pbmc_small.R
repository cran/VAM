### R code from vignette source 'VAM_pbmc_small.Rnw'

###################################################
### code chunk number 1: VAM_pbmc_small.Rnw:16-23
###################################################
library(VAM)
if (!requireNamespace("Seurat", quietly=TRUE)) {
   stop("Seurat package not available!")
}
if (!requireNamespace("SeuratObject", quietly=TRUE)) {
   stop("SeuratObject package not available!")
}


###################################################
### code chunk number 2: VAM_pbmc_small.Rnw:30-34
###################################################
SeuratObject::pbmc_small
gene.names = rownames(SeuratObject::pbmc_small)
gene.names[1:5]
Seurat::VariableFeatures(SeuratObject::pbmc_small)[1:5]


###################################################
### code chunk number 3: VAM_pbmc_small.Rnw:41-53
###################################################
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


###################################################
### code chunk number 4: VAM_pbmc_small.Rnw:60-63
###################################################
pbmc.vam = vamForSeurat(seurat.data=SeuratObject::pbmc_small,
    gene.set.collection=gene.set.collection,
    center=F, gamma=T, sample.cov=F, return.dist=T)


###################################################
### code chunk number 5: VAM_pbmc_small.Rnw:68-70
###################################################
pbmc.vam@assays$VAMdist[1,1:10]
pbmc.vam@assays$VAMcdf[1,1:10]


###################################################
### code chunk number 6: VAM_pbmc_small.Rnw:77-79
###################################################
Seurat::DefaultAssay(object = pbmc.vam) = "VAMcdf"
Seurat::FeaturePlot(pbmc.vam, reduction="tsne", features=gene.set.name)


