#' Get transcript annotation file
#'
#' @param organism either "mouse", "fly", "planaria", "worm", or "fish"
#' @export
#' @rdname mbl_get_gene_annotation
#' @importFrom readr read_tsv cols
#' @return transcript level annotations
mbl_get_transcript_annotation <- function(
  organism = c("mouse", "fly", "fish", "planaria", "worm", "mouse_mCherry")) {
  organism <- match.arg(organism)

  url <-paste0("https://s3.amazonaws.com/mbl.data/references/%s/",
               "gene_annotations.txt.gz")
  ti <- readr::read_tsv(sprintf(url, organism), col_types = readr::cols())
  colnames(ti) <- c("gene_id", "transcript_id", "gene_type", "symbol")[
    seq.int(ncol(ti))
  ]
  # rownames(ti) <- ti$target_id
  return(ti)
}

