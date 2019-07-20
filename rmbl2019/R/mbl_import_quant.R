#' Import salmon quantification results
#'
#' @export
#' @param path path to folder with quantification results
#' @param organism either "mouse", "fly", "planaria", "worm", or "fish"
#' @return a DGElist object
#' @importFrom checkmate assert_directory
#' @importFrom edgeR DGEList
#' @importFrom stringr str_match
#' @importFrom tximport tximport
mbl_import_quantitation_results <- function(
  path = path.expand("~/analysis/quantitation"),
  organism = c("mouse", "fly", "fish", "planaria", "worm"),
  dataset = c("mbl", "pre_mbl"),
  rm.description = TRUE) {
  # Ensure user asked for a valid organism
  organism <- match.arg(organism)
  dataset <- match.arg(dataset)
  checkmate::assert_directory(path)

  file_pattern <- switch(
    dataset,
    pre_mbl = ".*/(.*)_S\\d*_.*/quant.sf",
    mbl = ".*/(.*)/quant.sf"
  )
  files <- dir(path, pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
  names(files) <- stringr::str_match(
    string = files,
    pattern = file_pattern)[, 2]

  gene_anno <- mbl_get_transcript_annotation(organism)
  tx2gene <- gene_anno[, c("transcript_id", "gene_id")]

  txi <- tximport(
    files, type = "salmon",
    tx2gene = as.data.frame(gene_anno[, c("transcript_id", "gene_id")]),
    countsFromAbundance = "lengthScaledTPM")

  genes <- gene_anno[match(row.names(txi$counts), gene_anno$gene_id),
                     intersect(c("gene_id", "gene_type", "symbol"),
                               colnames(gene_anno))]

  # combine gene annotations & counts into a DGEList object
  y <- DGEList(txi$counts, genes = genes)
  return(y)
}
