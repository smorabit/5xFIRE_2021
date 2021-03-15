

library(Seurat)
library(DoubletFinder)


# Note: this should be done for each sample separately
meta <- data.frame()
for(cur_sample in unique(NucSeq$SampleID)){

  seurat_obj <- subset(NucSeq, SampleID == cur_sample)

  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list <- paramSweep_v3(seurat_obj, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)

  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  homotypic.prop <- modelHomotypic(seurat_obj$seurat_clusters)
  nExp_poi <- round(0.075*nrow(seurat_obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  seurat_obj <- doubletFinder_v3(seurat_obj, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  seurat_obj$Doublet <- seurat_obj@meta.data[[paste0('DF.classifications_0.25_0.09_', nExp_poi)]]

  # plot doublets on UMAP:
  # p <- DimPlot(seurat_obj,  group.by='Doublet', split.by='Doublet')
  # pdf(paste0(cur_sample, '_umap_doublets.pdf'), width=10, height=5)
  # p
  # dev.off()

  # update meta table:
  meta <- rbind(meta, seurat_obj@meta.data)

}
