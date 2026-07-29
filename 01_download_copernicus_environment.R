library(dplyr)
library(readr)
library(reticulate)
library(ncdf4)
library(purrr)


cm <- import("copernicusmarine")

ps <- readRDS("data/phyloseq_raw_clean.rds")
df <- data.frame(phyloseq::sample_data(ps), stringsAsFactors = FALSE)


df <- df |>
  transmute(
    sample_id,
    longitude,
    latitude,
    sampling_date
  )

datasets <- list(
  physics = list(
    id = "cmems_mod_glo_phy_myint_0.083deg_P1D-m",
    vars = c("thetao", "so", "uo", "vo", "zos"),
    depth = 0.49402499198913574
  ),
  bgc_model = list(
    id = "cmems_mod_glo_bgc_my_0.25deg_P1D-m",
    vars = c("chl", "o2", "no3", "po4", "si"),
    depth = 0.5057600140571594
  )
)

extract_nc_values <- function(ncfile, vars) {
  nc <- nc_open(ncfile)
  out <- list()
  for (v in vars) {
    out[[v]] <- if (v %in% names(nc$var)) as.numeric(ncvar_get(nc, v)) else NA
  }
  nc_close(nc)
  out
}

get_sample_data_single <- function(sample, ds) {
  tmpfile <- tempfile(fileext = ".nc")
  out <- tryCatch({
    cm$subset(
      dataset_id = ds$id,
      variables = ds$vars,
      minimum_longitude = as.numeric(sample[["longitude"]]),
      maximum_longitude = as.numeric(sample[["longitude"]]),
      minimum_latitude = as.numeric(sample[["latitude"]]),
      maximum_latitude = as.numeric(sample[["latitude"]]),
      start_datetime = format(as.Date(sample[["sampling_date"]]), "%Y-%m-%d"),
      end_datetime = format(as.Date(sample[["sampling_date"]]), "%Y-%m-%d"),
      minimum_depth = ds$depth,
      maximum_depth = ds$depth,
      coordinates_selection_method = "nearest",
      output_filename = tmpfile
    )
    extract_nc_values(tmpfile, ds$vars)
  }, error = function(e) {
    setNames(as.list(rep(NA, length(ds$vars))), ds$vars)
  })
  if (file.exists(tmpfile)) file.remove(tmpfile)
  out
}

get_sample_data <- function(sample) {
  result <- as.list(sample)
  for (ds_name in names(datasets)) {
    result <- c(result, get_sample_data_single(sample, datasets[[ds_name]]))
  }
  as.data.frame(result)
}

copernicus_env <- map_dfr(seq_len(nrow(df)), function(i) get_sample_data(df[i, ]))

