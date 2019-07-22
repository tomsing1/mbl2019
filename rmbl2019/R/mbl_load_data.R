#' Load DGEList with expression data from aws.s3
#'
#' @export
#' @param organism either "mouse", "fly", "planaria", "worm", or "fish"
#' @param dataset either "prem_mbl" or "mbl"
#' @importFrom janitor remove_empty
#' @return a DGElist object
mbl_load_data <- function(
  organism = c("mouse", "fly", "fish", "planaria", "worm"),
  dataset = c("pre_mbl", "mbl")) {
  organism <- match.arg(organism)
  dataset <- match.arg(dataset)
  if (organism == "fish" & dataset == "pre_mbl") {
    stop(paste(
      'Planarian data was not generated before the course at MBL!',
      'Maybe you need to use dataset = "mbl"?'))
  }
  if (organism == "planaria" & dataset == "mbl") {
    stop(paste(
      'Planarian data was not generated during the course at MBL!',
      'Maybe you need to use dataset = "pre_mbl"?'))
  }
  if (organism == "worm" & dataset == "mbl") {
    stop(paste(
      'C.elegans data was not generated during the course at MBL!',
      'Maybe you need to use dataset = "pre_mbl"?'))
  }
  s3_url <- sprintf("https://s3.amazonaws.com/mbl.data/dgelists/%s/%s.rds",
                 dataset, organism)
  y <- readRDS(url(s3_url))
  y$samples <- janitor::remove_empty(y$samples, which="cols")
  return(y)
}
