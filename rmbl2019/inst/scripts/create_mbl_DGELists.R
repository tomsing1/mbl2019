# Creating DGELists from libraries generated at MBL

library(purrr)
library(rmbl2019)
library(org.Mm.eg.db)  # mouse
library(org.Ce.eg.db)  # worm
library(org.Dr.eg.db)  # zebrafish
library(org.Dm.eg.db)  # fly
data("samples_pre_mbl")
row.names(samples_pre_mbl) <- samples_pre_mbl$sample_name

kUrl <- "s3://mbl.data/quantitation/pre_mbl/%s"
kAnnotations <- list(
  mouse = "org.Mm.eg.db",
  worm = "org.Ce.eg.db",
  fish = "org.Dr.eg.db",
  fly = "org.Dm.eg.db"
)
species_list <- c("mouse", "worm", "planaria", "fly")

add_entrez_ids <- function(dge, species, annotations = kAnnotations) {
  if (species %in% names(annotations)) {
    dge$genes$entrezid <- mapIds(
      get(annotations[[species]]), keys = row.names(dge), column = "ENTREZID",
      keytype = "ENSEMBL")
  } else {
    dge$genes$entrezid <- NA_character_
  }
  return(dge)
}

create_DGElist <- function(species) {
  message("Processing ", species, " data...")
  # download quantitation results from AWS S3
  temp_dir <-tempfile()
  dir.create(temp_dir, showWarnings = FALSE)
  aws.command <- paste("aws s3 sync", sprintf(kUrl, species), temp_dir)
  system(aws.command)

  # Create DGElist with gene annotations, but no sample annotations
  y <- mbl_import_quantitation_results(temp_dir, organism = species)

  # Add sample annotations
  not_annotated <- setdiff(colnames(y), samples_pre_mbl$sample_name)
  if (length(not_annotated) > 0) {
    message(sprintf("%s %s samples are not annotated in the sample sheet: %s",
                    length(not_annotated), species, paste(not_annotated, collapse = ", ")))
  }
  y <- y[, intersect(colnames(y), samples_pre_mbl$sample_name)]
  y$samples <- data.frame(
    y$samples,
    samples_pre_mbl[colnames(y), ]
  )
  row.names(y$samples) <- colnames(y)
  y <- add_entrez_ids(y, species = species, annotations = kAnnotations)
  unlink(temp_dir)
  return(y)
}

Main <- function() {
  # create all DGELists
  dges <- map(setNames(species_list, species_list), create_DGElist)

  # save DGELists as RDS files to the local hoome directory
  for (species in names(dges)) {
    saveRDS(dges[[species]], sprintf("~/%s.rds", species))
  }
  return(dges)
}

if (!interactive()) {
  dges <- Main()
}
