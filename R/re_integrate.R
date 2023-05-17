
#' Re Integrate
#'
#' @import
#' @param query
#' @param reference
#'
#' @return
#' @export
label_transfer <- function(query, reference) {





  ## **Subset, Split and re-integrate the atlas object**
  ## Subset the atlas data to only include celltypes involved in the Tcell exhaustion trajectory.
  # Subsetting the Atlas data into 3 types of Tcells for the exhaustion trajectory.
  atlas.seu <- subset(x = atlas.seu, subset = lv1_annot %in% c("T cells naive",
                                                               "CD8 cytotoxic",
                                                               "CD8 terminally exhausted"))


  ## Split the atlast object by source.
  atlas_sources <- SplitObject(atlas.seu, split.by = "source")

  ## Make an empty list to store the filtered seurat objects
  filtered_atlas_sources <- list()

  ## Loop through the `atlas_sources` and add objects which have more than 100 cells to the new list.
  for (i in 1:length(atlas_sources)) {
    # Check the number of cells in the Seurat object
    num_cells <- dim(atlas_sources[[i]]@assays$RNA@counts)[2]
    # If the number of cells is less than 100, skip this Seurat object
    if (num_cells < 100) {
      next
    }
    # Otherwise, add the Seurat object to the filtered list
    filtered_atlas_sources[[length(filtered_atlas_sources)+1]] <- atlas_sources[[i]]
  }

  ## Run SCTransform on each object in the list
  for (i in 1:length(filtered_atlas_sources)) {
    filtered_atlas_sources[[i]] <- SCTransform(filtered_atlas_sources[[i]], method = "glmGamPoi", verbose = FALSE)
  }

  #
  # Next, select features for downstream integration, and run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated.
  #
  # pancreas.features <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)
  # pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = pancreas.features,
  #                                     verbose = FALSE)
  # Next, identify anchors and integrate the datasets. Commands are identical to the standard workflow, but make sure to set normalization.method = 'SCT':
  #
  #   pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT",
  #                                              anchor.features = pancreas.features, verbose = FALSE)
  # pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT",
  #                                      verbose = FALSE)
  # Now proceed with downstream analysis (i.e. visualization, clustering) on the integrated dataset. Commands are identical to the standard workflow, but do not run the ScaleData function after integration. You can see that after integration, cells group by their biological cell type (which has been pre-annotated), instead of by their underlying technology.
  #
  # pancreas.integrated <- RunPCA(pancreas.integrated, verbose = FALSE)
  # pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30)
  # plots <- DimPlot(pancreas.integrated, group.by = c("tech", "celltype"), combine = FALSE)
  # plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3,
  #                                                                                                               byrow = TRUE, override.aes = list(size = 3))))
  # CombinePlots(plots)


  ## Intregrate all genes common between snRNA and scRNA.
  seuAnch <- FindIntegrationAnchors(filtered_atlas_sources, anchor.features = 1500,
                                    normalization.method = "LogNormalize",
                                    reduction = "cca", dims = 1:30,
                                    k.anchor = 5, k.filter = 200,
                                    k.score = 30, max.features = 200)

  atlas.seu <- IntegrateData(anchorset = seuAnch, dims = 1:30,
                             normalization.method = "LogNormalize")

  DefaultAssay(atlas.seu) = "integrated"
  #Should this be another SCTRANSFORM?
  #atlas.seu <- ScaleData(atlas.seu, verbose = FALSE)
  #atlas.seu <- FindVariableFeatures(atlas.seu)
  atlas.seu <- SCTransform(atlas.seu, method = "glmGamPoi", verbose = FALSE)
  atlas.seu <- RunPCA(atlas.seu, npcs = 30, verbose = FALSE)
  atlas.seu <- RunUMAP(atlas.seu, reduction = "pca", dims = 1:30)
  atlas.seu <- FindNeighbors(atlas.seu, reduction = "pca", dims = 1:30)
  atlas.seu <- FindClusters(atlas.seu, resolution = 0.5)


  ### Plot
  p1 <- DimPlot(atlas.seu, reduction = "umap", group.by = "lv1_annot")
  p2 <- DimPlot(atlas.seu, reduction = "umap", label = TRUE, repel = TRUE)
  p1 + p2






}



