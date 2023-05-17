
#' Title
#'
#' @param arg1
#' @param arg2
#'
#' @return
#' @export
#'
#' @examples
label_transfer <- function(arg1, arg2) {
  ## Find transfer anchors
  pre_anchors <- FindTransferAnchors(
    reference = atlas.seu,
    query = pre.seu,
    normalization.method = "SCT",
    reference.reduction = "pca",
    dims = 1:20
  )
  ## Map Query
  pre.seu <- MapQuery(
    anchorset = pre_anchors,
    query = pre.seu,
    reference = atlas.seu,
    refdata = list(
      cell_type = "lv1_annot"
    ),
    reference.reduction = "pca"
  )

  ## Remove cells from the patient data which don't have a confident cell_type prediction.
  pre.seu <- subset(pre.seu, subset = predicted.cell_type.score > 0.65)

  ## Process the patient data using Seurat
  #SCTransform again here?
  pre.seu <- RunPCA(pre.seu)
  pre.seu <- FindNeighbors(pre.seu)
  pre.seu <- FindClusters(pre.seu)
  pre.seu <- RunUMAP(pre.seu, dims = 1:30, n_neighbors = 20, min_dist = 0.3)
}


