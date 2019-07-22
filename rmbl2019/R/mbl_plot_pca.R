#' Plot the first principal compontents for a DGEList
#'
#' @param x DGEList
#' @param intgroup Character vector with column names of x$samples to color
#' the points
#' @param prior.count Number of pseudocounts to add before log transformation
#' @param ntop Number of most variable genes to use for PCA
#' @return ggplot2 object
#' @importFrom matrixStats rowVars
#' @export
#' @importFrom ggplot2 ggplot geom_point xlab ylab coord_fixed aes_string
#' @note This function was modified from the DESeq2::plotPCA function
mbl_plot_pca <- function(x, intgroup = NA, prior.count = 5, ntop = 500) {
  cpms <- cpm(x, log = TRUE, prior.count = prior.count)
  rv <- matrixStats::rowVars(cpms)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(cpms[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!is.na(intgroup)) {
    if (!all(intgroup %in% names(x$samples))) {
      stop("the argument 'intgroup' should specify columns of x$samples")
    }
    intgroup.df <- as.data.frame(x$samples[, intgroup, drop = FALSE])
    if (length(intgroup) > 1) {
      group <- factor(apply(intgroup.df, 1, paste, collapse = ":"))
    } else {
      group <- intgroup.df[[intgroup]]
    }
    d <- data.frame(
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      group,
      name = colnames(x))
    p <- ggplot2::ggplot(data = d,
                         ggplot2::aes_string(x = "PC1", y = "PC2",
                                            color = "group")) +
      ggplot2::geom_point(size = 3) +
      ggplot2::xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
      ggplot2::ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
      ggplot2::coord_fixed()
  } else {
    d <- data.frame(
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      name = colnames(x))
    p <- ggplot2::ggplot(data = d, aes_string(x = "PC1", y = "PC2")) +
      ggplot2::geom_point(size = 3) +
      ggplot2::xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
      ggplot2::ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
      ggplot2::coord_fixed()
  }
    return(p)
}
