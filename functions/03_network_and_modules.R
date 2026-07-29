library(igraph)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(vegan)
library(RColorBrewer)

module_abundance <- function(ps, relative = TRUE) {
  cl <- readRDS("output/raw_tables/cluster.rds") |>
    transmute(ASV = OTU_ID, Module = as.character(Cluster))
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  cl <- cl |> filter(ASV %in% rownames(otu))
  otu <- otu[cl$ASV, , drop = FALSE]
  m <- rowsum(otu, group = cl$Module)
  if (relative) m <- sweep(m, 2, colSums(m), "/")
  ord <- suppressWarnings(order(as.numeric(rownames(m)), rownames(m), na.last = TRUE))
  m[ord, , drop = FALSE]
}

top_modules <- function(module_mat, n = 4) {
  names(sort(rowMeans(module_mat, na.rm = TRUE), decreasing = TRUE))[seq_len(n)]
}

module_correlations <- function(module_mat, meta) {
  env <- c("sst_degC", "salinity_psu", "chlorophyll_mg_m3", "oxygen_mmol_m3", "nitrate_mmol_m3", "phosphate_mmol_m3", "silicate_mmol_m3", "longitude", "latitude")
  long <- as.data.frame(t(module_mat)) |> rownames_to_column("sample_id")
  dat <- left_join(long, meta, by = "sample_id")
  bind_rows(lapply(rownames(module_mat), function(m) {
    bind_rows(lapply(env, function(v) {
      z <- suppressWarnings(cor.test(dat[[m]], as.numeric(dat[[v]]), method = "spearman", exact = FALSE))
      data.frame(Module = m, Variable = v, rho = unname(z$estimate), p = z$p.value)
    }))
  })) |>
    mutate(FDR = p.adjust(p, method = "BH"))
}

weighted_sd <- function(x, w) {
  m <- weighted.mean(x, w, na.rm = TRUE)
  sqrt(weighted.mean((x - m)^2, w, na.rm = TRUE))
}

module_thermal_niche <- function(module_mat, meta, modules) {
  temp <- as.numeric(meta$sst_degC[match(colnames(module_mat), meta$sample_id)])
  bind_rows(lapply(modules, function(m) {
    w <- as.numeric(module_mat[m, ])
    data.frame(Module = m, optimum = weighted.mean(temp, w, na.rm = TRUE), breadth = weighted_sd(temp, w), mean_abundance = mean(w, na.rm = TRUE))
  }))
}

network_layout <- function(g) {
  set.seed(123456789)
  x <- layout_with_fr(g, niter = 5000, grid = "nogrid")
  x <- scale(x)
  lon_name <- if ("longitude" %in% vertex_attr_names(g)) "longitude" else "Longitude.x"
  lon <- as.numeric(vertex_attr(g, lon_name))
  if (cor(x[, 1], lon, use = "complete.obs") < 0) x[, 1] <- -x[, 1]
  x * 0.03
}

network_node_size <- function(g) {
  d <- degree(g, mode = "all")
  ifelse(log(pmax(d, 1)) < 3, 2.5, log(pmax(d, 1)) * 1.8)
}

network_continuous_colors <- function(x, palette, reverse = FALSE) {
  pal <- brewer.pal(9, palette)
  if (reverse) pal <- rev(pal)
  z <- vegan::decostand(as.numeric(x), method = "range")
  rgb(colorRamp(pal)(z), maxColorValue = 255)
}

network_continuous_scale <- function(x, palette, reverse = FALSE, label = NULL) {
  x <- as.numeric(x)
  rng <- range(x, na.rm = TRUE)
  pal <- brewer.pal(9, palette)
  if (reverse) pal <- rev(pal)
  cols <- colorRampPalette(pal)(100)
  usr <- par("usr")
  dx <- diff(usr[1:2])
  dy <- diff(usr[3:4])
  x0 <- usr[2] - 0.07 * dx
  x1 <- usr[2] - 0.035 * dx
  y0 <- usr[3] + 0.10 * dy
  y1 <- usr[3] + 0.46 * dy
  ys <- seq(y0, y1, length.out = 101)
  rect(x0, ys[-101], x1, ys[-1], col = cols, border = cols, xpd = NA)
  rect(x0, y0, x1, y1, border = "grey30", lwd = 0.8, xpd = NA)
  text(x1 + 0.01 * dx, c(y0, y1), round(rng, 2), adj = c(0, 0.5), cex = 0.65, xpd = NA)
  if (!is.null(label)) text((x0 + x1) / 2, y1 + 0.05 * dy, label, cex = 0.72, xpd = NA)
}
