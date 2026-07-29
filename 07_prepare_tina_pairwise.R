## TINA pairwise distance
# grazie dermilke
source("functions/01_atlantic_transect_and_community.R")

library(phyloseq)
library(vegan)
library(geosphere)
library(readr)
library(dplyr)
library(tibble)
library(ape)

tina_functions <- "data/tina/Functions_Similarity_Indices.R"


source(tina_functions)

ps <- read_atlantic_ps("css")
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
meta <- sample_metadata(ps)
sample_ids <- sample_names(ps)

tina <- distance_wrapper2(
  ps,
  method = "TINA_weighted",
  size.thresh = 1,
  pseudocount = 1e-6,
  nblocks = 25,
  use.cores = 2,
  cor.use = "na.or.complete"
)

rownames(tina) <- colnames(tina) <- sample_ids

bray <- as.matrix(distance(ps, method = "bray"))
coords <- as.matrix(meta[, c("longitude", "latitude")])
rownames(coords) <- sample_ids
geo <- distm(coords, fun = distHaversine) / 1000
rownames(geo) <- colnames(geo) <- sample_ids

lower_matrix <- function(mat, name) {
  idx <- which(lower.tri(mat), arr.ind = TRUE)
  out <- tibble(
    sample_1 = rownames(mat)[idx[, 1]],
    sample_2 = colnames(mat)[idx[, 2]],
    value = mat[idx]
  )
  names(out)[3] <- name
  out
}

tina_pairwise <- lower_matrix(as.matrix(as.dist(tina)), "tina_dist") |>
  left_join(lower_matrix(geo, "geo_km"), by = c("sample_1", "sample_2")) |>
  left_join(lower_matrix(bray, "bray"), by = c("sample_1", "sample_2"))

write_csv(tina_pairwise, "output/raw_tables/tina_pina_pairwise.csv")
tina_pairwise
