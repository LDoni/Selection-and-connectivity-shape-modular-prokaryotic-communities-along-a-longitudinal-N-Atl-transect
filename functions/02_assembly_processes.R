library(dplyr)
library(readr)
library(vegan)
library(geosphere)

classify_process <- function(bnti, rc) {
  case_when(
    bnti <= -2 ~ "Homogeneous Selection",
    bnti >= 2 ~ "Heterogeneous Selection",
    abs(bnti) < 2 & rc >= 0.95 ~ "Dispersal Limitation",
    abs(bnti) < 2 & rc <= -0.95 ~ "Homogenising Dispersal",
    TRUE ~ "Drift"
  )
}

read_assembly_pairs <- function() {
  read_csv("output/raw_tables/assembly_module_pairwise.csv", show_col_types = FALSE)
}
