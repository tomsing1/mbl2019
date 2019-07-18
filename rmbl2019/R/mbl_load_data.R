#' Load DGEList with expression data from aws.s3
#'
#' @export
#' @param organism either "mouse", "fly", "planaria", "worm", or "fish"
#' @param dataset either "prem_mbl" or "mbl"
#' @return a DGElist object
mbl_load_data <- function(
  organism = c("mouse", "fly", "fish", "planaria", "worm"),
  dataset = c("pre_mbl", "mbl")) {
  organism <- match.arg(organism)
  dataset <- match.arg(dataset)
  if (dataset == "mbl") {
    stop("Data for samples generated at MBL is not yet available, sorry.")
  }
  if (organism == "fish" & dataset == "pre_mbl") {
    stop("Fish data was not generated pre-MBL!")
  }
  s3_url <- sprintf("https://s3.amazonaws.com/mbl.data/dgelists/%s/%s.rds",
                 dataset, organism)
  y <- readRDS(url(s3_url))
  return(y)
}
