# This script copies and renames FASTQ files from a temporary S3 bucket into
# the final destination, e.g. subfolders for each species.

library(readr)
library(aws.s3)
library(here)
library(janitor)
library(dplyr)

kSampleAnno <- system.file("extdata", "mbl_samples.txt",
                           package = "rmbl2019")
kFileNames <- system.file("extdata", "mbl_seq_file_names.txt",
                          package = "rmbl2019")
kS3Bucket <- "mbl.data" # no trailing slash!
kS3InPrefix <- "reads/tmp"  # no trailing slash!
kS3OutPrefix <- "reads/mbl"  # no trailing slash!

Main <- function() {
  sample_anno <- read_tsv(kSampleAnno, col_types = cols()) %>%
    janitor::remove_empty(which = c("rows", "cols")) %>%
    janitor::clean_names()
  file_names <- read_tsv(kFileNames, col_types = cols()) %>%
    janitor::remove_empty(which = c("rows", "cols")) %>%
    janitor::clean_names()
  stopifnot(!any(duplicated(file_names$original_file_name)))
  stopifnot(!any(duplicated(file_names$target_file_name)))
  unannotated_fastq <- dplyr::anti_join(
    file_names, sample_anno, by = "sample_name")
  if (nrow(unannotated_fastq) > 0) {
    message(sprintf("%s FASTQ files were not annotated in %s",
                    nrow(unannotated_fastq), basename(kSampleAnno)))
    print(unannotated_fastq)
  }

  available_fastq <- aws.s3::get_bucket_df(bucket = kS3Bucket,
                                           prefix = kS3InPrefix) %>%
    dplyr::filter(grepl("fastq.gz$", Key)) %>%
    dplyr::mutate(Key = basename(Key)) %>%
    dplyr::select(Key)
  message(sprintf("%s source FASTQ files found.", nrow(available_fastq)))

  unexpected_fastq <- setdiff(
    available_fastq$Key, file_names$original_file_name)
  message(sprintf("%s FASTQ files were not listed in %s",
                  length(unexpected_fastq), basename(kFileNames)))

  purrr::walk(
    available_fastq$Key, .f = function(old_file) {
      message(sprintf("Processing existing FASTQ file %s", old_file))
      # identify new filename
      new_file_anno <- file_names %>%
        dplyr::filter(original_file_name == old_file)
      new_file <- new_file_anno$target_file_name

      # identify species
      new_sample_anno <- dplyr::inner_join(sample_anno, new_file_anno,
                                           by = "sample_name")
      if(!nrow(new_sample_anno) == 1L) {
        message("Skipping file %s because it is not annotated (with a species)",
                old_file)
        return(NULL)
      }
      species = switch(new_sample_anno$organism,
                       drosophila = "fly",
                       mouse = "mouse",
                       zebrafish = "fish",
                       NA_character_
      )
      if (is.na(species)) {
        stop(sprintf("Species %s wasn't recognized.", species))
      }
      message(sprintf("Copying file %s to %s", old_file, new_file))
      aws.s3::copy_object(from_object = paste0(kS3InPrefix, "/", old_file),
                          to_object = paste0(kS3OutPrefix, "/", species, "/",
                                             new_file),
                          from_bucket = kS3Bucket, to_bucket = kS3Bucket)
  })
}

if (!interactive()) {
  Main()
}
