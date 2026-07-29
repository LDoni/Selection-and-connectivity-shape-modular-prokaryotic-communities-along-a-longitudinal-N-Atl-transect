## Assembly pairwise table
#grazie dermilke 
source("functions/01_atlantic_transect_and_community.R")
source("functions/02_assembly_processes.R")
source("functions/03_network_and_modules.R")  

library(dplyr)
library(readr)
library(tibble)
library(vegan)

bnti_file <- "data/assembly/Prokaryotes_Atlantic_weighted_bNTI.csv"
rc_file <- "data/assembly/Raup_Crick_Prok.csv"

if (file.exists(bnti_file) && file.exists(rc_file)) {
  bnti <- read.csv(bnti_file, row.names = 1, check.names = FALSE)
  rc <- read.csv(rc_file, row.names = 1, check.names = FALSE)

  idx <- which(upper.tri(as.matrix(bnti)), arr.ind = TRUE)
  stegen <- data.frame(
    Sample_ID = rownames(bnti)[idx[, 1]],
    To_Sample = colnames(bnti)[idx[, 2]],
    bNTI = as.matrix(bnti)[idx],
    RC_BC = as.matrix(rc)[idx]
  ) |>
    mutate(Mechanism = classify_process(bNTI, RC_BC))
} else {
  stegen <- read_csv("output/raw_tables/stegen_pairwise.csv", show_col_types = FALSE)
}

ps <- read_atlantic_ps("css")
meta <- sample_metadata(ps) |>
  transmute(
    sample_id,
    Zone = as.character(oceanic_sector),
    Temp = sst_degC,
    Lon = longitude,
    Lat = latitude
  )

modules <- module_abundance(ps, relative = TRUE)
main_modules <- top_modules(modules, 4)
module_dist <- as.matrix(vegdist(t(modules[main_modules, , drop = FALSE]), method = "bray"))

assembly <- stegen |>
  rename(Sample_1 = Sample_ID, Sample_2 = To_Sample) |>
  rowwise() |>
  mutate(module_bray = module_dist[Sample_1, Sample_2]) |>
  ungroup() |>
  left_join(meta, by = c("Sample_1" = "sample_id")) |>
  rename(Zone_1 = Zone, Temp_1 = Temp, Lon_1 = Lon, Lat_1 = Lat) |>
  left_join(meta, by = c("Sample_2" = "sample_id")) |>
  rename(Zone_2 = Zone, Temp_2 = Temp, Lon_2 = Lon, Lat_2 = Lat) |>
  mutate(
    Pair_A = ifelse(Sample_1 <= Sample_2, Sample_1, Sample_2),
    Pair_B = ifelse(Sample_1 <= Sample_2, Sample_2, Sample_1),
    module_similarity = 1 - module_bray,
    temp_diff = abs(Temp_1 - Temp_2),
    lon_diff = abs(Lon_1 - Lon_2),
    geographic_km = geosphere::distHaversine(
      cbind(Lon_1, Lat_1), cbind(Lon_2, Lat_2)
    ) / 1000,
    zone_pair = ifelse(Zone_1 <= Zone_2, paste(Zone_1, Zone_2, sep = "-"), paste(Zone_2, Zone_1, sep = "-")),
    temp_bin = cut(temp_diff, breaks = c(0, 1, 2, 4, 6, 8), right = FALSE)
  ) |>
  select(
    Sample_1, Sample_2, bNTI, RC_BC, Mechanism, Pair_A, Pair_B,
    module_bray, Zone_1, Temp_1, Lon_1, Zone_2, Temp_2, Lon_2,
    module_similarity, temp_diff, lon_diff, geographic_km, zone_pair, temp_bin
  )

head(assembly)

write_csv(assembly, "output/raw_tables/assembly_module_pairwise.csv")

