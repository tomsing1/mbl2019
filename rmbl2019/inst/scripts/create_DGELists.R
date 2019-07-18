# Creating DGELists for pre-MBL datasets
# (This script should be executed on an AWS EC2 instance.)

library(purrr)
library(rmbl2019)
data("samples_pre_mbl")
row.names(samples_pre_mbl) <- samples_pre_mbl$sample_name

kUrl <- "s3://mbl.data/quantitation/pre_mbl/%s"
species_list <- c("mouse", "worm", "fly", "planaria")

create_DGElist <- function(species) {
  # download quantitation results from AWS S3
  temp_dir <-tempfile()
  dir.create(temp_dir, showWarnings = FALSE)
  aws.command <- paste("aws s3 sync", sprintf(kUrl, species), temp_dir)
  system(aws.command)

  # Create DGElist with gene annotations, but no sample annotations
  y <- mbl_import_quantitation_results(temp_dir, organism = species)

  # Add sample annotations
  stopifnot(all(colnames(y) %in% samples_pre_mbl$sample_name))
  y$samples <- data.frame(
    y$samples,
    samples_pre_mbl[colnames(y), ]
  )
  unlink(temp_dir)
  return(y)
}

# create all DGELists
dges <- map(setNames(species_list, species_list), create_DGElist)

# save DGELists as RDS files to the local hoome directory
for (species in names(dges)) {
  saveRDS(y, sprintf("~/%s.rds", species))
}
