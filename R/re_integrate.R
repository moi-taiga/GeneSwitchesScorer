
#' Re Integrate
#'
#' @import
#' @param seu
#'
#' @return
#' @export
re_integrate <- function(seu) {
  ## **Split and re-integrate the atlas object**

  ## Split the atlast object by source.
  split.seu <- SplitObject(seu, split.by = "source")

  ## Loop (backwards) through the `split.seu` and remove objects which have less than 100 cells.
  # Start the loop from the last index and iterate backward to the first index
  for (i in length(split.seu):1) {
    # Get the number of cells in the Seurat object
    num_cells <- dim(split.seu[[i]]@assays$RNA@counts)[2]

    # Check if the number of cells is less than 100
    if (num_cells < 100) {
      # Remove the Seurat object from the list
      split.seu <- split.seu[-i]
    }
  }

  ## Run SCTransform on each object in the list
  for (i in 1:length(split.seu)) {
    split.seu[[i]] <- SCTransform(split.seu[[i]], method = "glmGamPoi", verbose = FALSE)
  }

  ## Run SCTransform on each object in the list
  for (i in 1:length(atlas_sources)) {
    tryCatch({
      atlas_sources[[i]] <- SCTransform(atlas_sources[[i]], method = "glmGamPoi", verbose = FALSE)
    }, error = function(e) {
      message <- paste("Error occurred while applying SCTransform to Seurat object", i, ":")
      message <- paste0(message,
                       e$message,
                       "\n",
                       "To resolve this error Use \`fixInNamespace(\"make_cell_attr\", \"sctransform\")\` to change !identical to !setequal.",
                       "\n",
                       "This is expected to be default in the next update of Seurat",
                       "\n",
                       "\n")
      cat(message)
    })
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
  seuAnch <- FindIntegrationAnchors(filtered_split.seu, anchor.features = 1500,
                                    normalization.method = "LogNormalize",
                                    reduction = "cca", dims = 1:30,
                                    k.anchor = 5, k.filter = 200,
                                    k.score = 30, max.features = 200)

  seu <- IntegrateData(anchorset = seuAnch, dims = 1:30,
                             normalization.method = "LogNormalize")

  DefaultAssay(seu) = "integrated"
  #Should this be another SCTRANSFORM?
  #seu <- ScaleData(seu, verbose = FALSE)
  #seu <- FindVariableFeatures(seu)
  seu <- SCTransform(seu, method = "glmGamPoi", verbose = FALSE)
  seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:30)
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:30)
  seu <- FindClusters(seu, resolution = 0.5)


  ### Plot
  p1 <- DimPlot(seu, reduction = "umap", group.by = "lv1_annot")
  p2 <- DimPlot(seu, reduction = "umap", label = TRUE, repel = TRUE)
  p1 + p2






}



