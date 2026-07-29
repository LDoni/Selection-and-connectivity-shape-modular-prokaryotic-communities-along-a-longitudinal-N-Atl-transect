library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(geosphere)
library(scales)
library(metagenomeSeq)

zone_cols <- c(ANW = "#4280fc", ANC = "#ffb452", ANE = "#f7170a")

prepare_atlantic_ps <- function(ps) {
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu_table(ps) <- otu_table(t(otu), taxa_are_rows = TRUE)
  sample_data(ps)$oceanic_sector <- factor(as.character(sample_data(ps)$oceanic_sector), levels = names(zone_cols))
  ps
}

read_atlantic_raw_ps <- function() {
  prepare_atlantic_ps(readRDS("data/phyloseq_raw_clean.rds"))
}

css_normalize_ps <- function(ps) {
  ps <- prepare_atlantic_ps(ps)
  keep_taxa <- genefilter_sample(ps, filterfun_sample(function(x) x > 1), A = 1)
  ps <- prune_taxa(keep_taxa, ps)

  data_metaseq <- phyloseq_to_metagenomeSeq(ps)
  p <- metagenomeSeq::cumNormStat(data_metaseq)
  data_cumnorm <- metagenomeSeq::cumNorm(data_metaseq, p = p)
  data_css <- metagenomeSeq::MRcounts(data_cumnorm, norm = TRUE, log = TRUE)

  otu_table(ps) <- otu_table(as.matrix(data_css), taxa_are_rows = TRUE)
  prepare_atlantic_ps(ps)
}

tss_normalize_ps <- function(ps) {
  ps <- prepare_atlantic_ps(ps)
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) {
    otu <- t(otu)
  }
  otu <- sweep(otu, 2, colSums(otu), "/")
  otu_table(ps) <- otu_table(otu, taxa_are_rows = TRUE)
  prepare_atlantic_ps(ps)
}

read_atlantic_ps <- function(normalization = "css") {
  normalization <- match.arg(normalization, c("css", "raw", "tss"))
  ps <- read_atlantic_raw_ps()
  if (normalization == "raw") return(ps)
  if (normalization == "tss") return(tss_normalize_ps(ps))
  css_normalize_ps(ps)
}

sample_metadata <- function(ps) {
  data.frame(sample_data(ps), stringsAsFactors = FALSE) |>
    rownames_to_column("sample_id_from_rownames") |>
    mutate(
      oceanic_sector = factor(oceanic_sector, levels = names(zone_cols))
    )
}

geo_matrix <- function(meta) {
  xy <- as.matrix(meta[, c("longitude", "latitude")])
  out <- distm(xy, fun = distHaversine) / 1000
  rownames(out) <- colnames(out) <- meta$sample_id
  out
}

lower_df <- function(mat, value) {
  idx <- which(lower.tri(mat), arr.ind = TRUE)
  data.frame(
    Sample_1 = rownames(mat)[idx[, 1]],
    Sample_2 = colnames(mat)[idx[, 2]],
    value = mat[idx],
    stringsAsFactors = FALSE
  ) |>
    rename(!!value := value)
}
