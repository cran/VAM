### R code from vignette source 'VAM_pbmc_small_sctransform.Rnw'

###################################################
### code chunk number 1: VAM_pbmc_small_sctransform.Rnw:15-19
###################################################
library(VAM)
if (!requireNamespace("Seurat", quietly=TRUE)) {
   stop("Seurat package not available!")
}


###################################################
### code chunk number 2: VAM_pbmc_small_sctransform.Rnw:26-29
###################################################
SeuratObject::pbmc_small
gene.names = rownames(SeuratObject::pbmc_small)
gene.names[1:5]


###################################################
### code chunk number 3: VAM_pbmc_small_sctransform.Rnw:34-39
###################################################
pbmc_sctransform = Seurat::SCTransform(SeuratObject::pbmc_small, verbose=F)
# Compute PCA and UMAP on the normalized values
pbmc_sctransform = Seurat::RunPCA(pbmc_sctransform, npcs=10)
pbmc_sctransform = Seurat::RunUMAP(pbmc_sctransform, dims = 1:10) 
Seurat::VariableFeatures(pbmc_sctransform)[1:5]


###################################################
### code chunk number 4: VAM_pbmc_small_sctransform.Rnw:46-58
###################################################
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


###################################################
### code chunk number 5: VAM_pbmc_small_sctransform.Rnw:65-68
###################################################
pbmc.vam = vamForSeurat(seurat.data=pbmc_sctransform,
    gene.set.collection=gene.set.collection,
    center=F, gamma=T, sample.cov=F, return.dist=T)


###################################################
### code chunk number 6: VAM_pbmc_small_sctransform.Rnw:73-75
###################################################
pbmc.vam@assays$VAMdist[1,1:10]
pbmc.vam@assays$VAMcdf[1,1:10]


###################################################
### code chunk number 7: VAM_pbmc_small_sctransform.Rnw:82-84
###################################################
Seurat::DefaultAssay(object = pbmc.vam) = "VAMcdf"
Seurat::FeaturePlot(pbmc.vam, reduction="tsne", features=gene.set.name)


