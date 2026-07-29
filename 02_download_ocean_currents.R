library(reticulate)
library(dplyr)

#data for oceanographic models

cm <- import("copernicusmarine")

ps <- readRDS("data/phyloseq_raw_clean.rds")
stations <- data.frame(phyloseq::sample_data(ps), stringsAsFactors = FALSE)

dir.create("data", showWarnings = FALSE)

stations |>
  transmute(
    station = sample_id,
    date = sampling_date,
    lat = latitude,
    lon = longitude
  ) |>
  write.csv("data/stations.csv", row.names = FALSE)

cm$subset(
  dataset_id = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
  variables = c("uo", "vo"),
  minimum_longitude = -100,
  maximum_longitude = 20,
  minimum_latitude = 0,
  maximum_latitude = 70,
  start_datetime = "2020-11-23",
  end_datetime = "2023-12-26",
  minimum_depth = 0,
  maximum_depth = 1,
  output_filename = "data/ocean_currents_surface_2020_2023.nc"
)
