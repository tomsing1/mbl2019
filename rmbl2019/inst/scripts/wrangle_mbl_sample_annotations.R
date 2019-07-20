library(readr)
library(rmbl2019)
library(janitor)
library(dplyr)
library(usethis)

samples_mbl <- readr::read_tsv(
  file = system.file("extdata", "mbl_samples.txt", package = "rmbl2019"),
  col_types = readr::cols()) %>%
  janitor::remove_empty(which = c("rows", "cols")) %>%
  janitor::clean_names() %>%
  dplyr::mutate(organism = case_when(
    organism == "drosophila" ~ "fly",
    organism == "mouse" ~ "mouse",
    organism == "zebrafish" ~ "fish",
    TRUE ~ NA_character_
  ))

usethis::use_data(samples_mbl, overwrite = TRUE)
