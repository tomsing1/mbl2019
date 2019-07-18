#' Create tidy table with expression data from DGEList
#'
#' This function returns the log2(CPM + 1) values for all genes and samples
#' in the form of a tidy data.frame.
#' @export
#' @param y DGElist
#' @param log_transform Log tranform the CPMs?
#' @param prior.count Prior count to add to the CPMs
#' @return a tidy data.frame
#' @importFrom tibble rownames_to_column
#' @importFrom edgeR cpm
#' @importFrom dplyr inner_join everything select
#' @importFrom tidyr gather
mbl_create_tidy_table <- function(y, log_transform = TRUE, prior.count = 1) {
  if (!is(y, "DGEList")) {
    stop("Object y must be a DGEList object!")
  }
  fdata <- y$genes
  pdata <- tibble::rownames_to_column(y$samples, "sample_id")
  if (log_transform == TRUE) {
    message(
      sprintf("Returning log2 transformed CPMs with %s pseudocount(s)!",
              prior.count))
  }
  cpms <- edgeR::cpm(y, log = log_transform,
                     prior.count = prior.count) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    tidyr::gather(key = "sample_id", value = "cpm", -gene_id) %>%
    dplyr::inner_join(pdata, by = "sample_id") %>%
    dplyr::inner_join(fdata, by = "gene_id") %>%
    dplyr::select(symbol, gene_id, sample_id, everything())
  return(cpms)
}
