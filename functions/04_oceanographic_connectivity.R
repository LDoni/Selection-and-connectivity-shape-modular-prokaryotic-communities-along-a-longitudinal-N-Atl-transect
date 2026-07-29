library(readr)
library(dplyr)

read_ocean_pairs <- function() {
  read_csv("output/raw_tables/oceanparcels_stationdates_bidirectional_pairwise_365d50km7.csv", show_col_types = FALSE)
}

read_forward_connectivity <- function() {
  read.csv("data/connectivity/forward_connectivity_matrix.csv", row.names = 1, check.names = FALSE)
}

read_backward_connectivity <- function() {
  read.csv("data/connectivity/backward_connectivity_matrix.csv", row.names = 1, check.names = FALSE)
}

read_forward_travel_time <- function() {
  read.csv("data/connectivity/forward_min_travel_days.csv", row.names = 1, check.names = FALSE)
}

read_backward_travel_time <- function() {
  read.csv("data/connectivity/backward_min_travel_days.csv", row.names = 1, check.names = FALSE)
}
