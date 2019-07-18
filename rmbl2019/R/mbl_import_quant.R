#' Import salmon quantification results
#'
#' @export
#' @param path path to folder with quantification results
#' @param organism either "mouse", "fly", "planaria", "worm", or "fish"
#' @return a DGElist object
#' @importFrom checkmate assert_directory
#' @importFrom edgeR DGEList
#' @importFrom stringr str_match
mbl_import_quantitation_results <- function(
  path,
  organism = c("mouse", "fly", "fish", "planaria", "worm"),
  rm.description = TRUE) {
  # Ensure user asked for a valid organism
  organism <- match.arg(organism)
  checkmate::assert_directory(path)

  files <- dir(path, pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
  names(files) <- str_match(files, pattern = ".*/(.*)_S\\d*_.*/quant.sf")[, 2]

  gene_anno <- mbl_get_transcript_annotation(organism)
  tx2gene <- gene_anno[, c("transcript_id", "gene_id")]

  txi <- tximport(
    files, type = "salmon",
    tx2gene = as.data.frame(gene_anno[, c("transcript_id", "gene_id")]),
    countsFromAbundance = "lengthScaledTPM")

  genes <- gene_anno[match(row.names(txi$counts), gene_anno$gene_id),
                     c("gene_id", "gene_type", "symbol")]

  # combine gene annotations & counts into a DGEList object
  y <- DGEList(txi$counts, genes = genes)
  return(y)
}
