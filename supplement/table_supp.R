library(tidyverse)

#Combine phenotype data and growth data for all pure isolates.

GRBC_pure_isolate_taxonomy <- read_csv("fig1/GRBC_pure_isolate_taxonomy.csv")

GRBC_collection_phenotype <- read_csv("Supplement/GRBC_collection_phenotype.csv", 
                                      col_types = cols(`D3/D31 Number` = col_character(), 
                                                       `GRBC Collection Number` = col_character())) |>
  rename("capture_day" = ...1)
