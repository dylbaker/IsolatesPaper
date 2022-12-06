library(tidyverse)
library(Biostrings)
library(gridExtra)

#Combine phenotype data and growth data for all pure isolates.

GRBC_pure_isolate_taxonomy <- read_csv("fig1/GRBC_pure_isolate_taxonomy.csv") |>
  mutate(Group = str_replace(Group, "point", "p"),
         Isolate_Number = str_replace(Isolate_Number, "_D31", "D31"),
         Isolate_Number = case_when(Isolate_Number == "431_D3" ~ "43_1D3",
                                    Isolate_Number == "432_D3" ~ "43_2D3",
                                    T ~ Isolate_Number)) |>
  select(-...1)

GRBC_collection_phenotype <-
  read_csv(
    "Supplement/GRBC_collection_phenotype.csv",
    col_types = cols(
      `D3/D31 Number` = col_character(),
      `GRBC Collection Number` = col_character())
  ) |>
  dplyr::rename(
    "time" = ...1,
    "number" = `D3/D31 Number`,
    "orig_morph" = `Original Morphology`,
    "updated_morph" = `Most Recent Morphology (on plate)`,
    "freezer_morph" = `Morphology of Freezer Stock Restreak`
  ) |>
  dplyr::slice(1:463) |>
  mutate(number = str_replace_all(number, "\\.", "_")) |>
  unite(col = "Isolate_Number", number:time, sep = '') |>
  select(Isolate_Number,
         `Algal Host`,
         orig_morph,
         updated_morph,
         freezer_morph) |>
  #Fill in missing morph data
  mutate(
    updated_morph = case_when(is.na(updated_morph) ~ freezer_morph,
                              T ~ updated_morph),
    updated_morph = case_when(is.na(updated_morph) ~ orig_morph,
                              T ~ updated_morph)
  )
  
fasta <- readDNAStringSet("mothur_outputs/all_seqs/asv_list.fasta")
asv <- names(fasta)
rep_seq <- paste(fasta)
rep_seqs <- data.frame(asv, rep_seq) |>
  mutate(asv = as.numeric(str_replace(asv, "ASV0+|ASV", "")))

GRBC_all_data <- left_join(GRBC_pure_isolate_taxonomy, GRBC_collection_phenotype) |>
  left_join(rep_seqs, by = "asv") |>
  select(-host_species)|>
  dplyr::rename("capture_day" = "time")

write.csv(GRBC_all_data, "Supplement/table_S1.csv")

