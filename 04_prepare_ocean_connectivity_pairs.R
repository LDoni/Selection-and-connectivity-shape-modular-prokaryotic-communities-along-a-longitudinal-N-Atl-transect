library(readr)
library(dplyr)
library(tidyr)

# for oceanographic analysis

forward_connectivity <- read.csv(
  "data/connectivity/forward_connectivity_matrix.csv",
  row.names = 1, check.names = FALSE
)
backward_connectivity <- read.csv(
  "data/connectivity/backward_connectivity_matrix.csv",
  row.names = 1, check.names = FALSE
)
forward_time <- read.csv(
  "data/connectivity/forward_min_travel_days.csv",
  row.names = 1, check.names = FALSE
)
backward_time <- read.csv(
  "data/connectivity/backward_min_travel_days.csv",
  row.names = 1, check.names = FALSE
)

stations <- rownames(forward_connectivity)
station_pairs <- combn(stations, 2, simplify = FALSE)

ocean_pairs <- lapply(station_pairs, function(x) {
  a <- x[1]
  b <- x[2]
  connectivity <- c(
    backward_connectivity[a, b], backward_connectivity[b, a],
    forward_connectivity[a, b], forward_connectivity[b, a]
  )
  travel_time <- c(
    backward_time[a, b], backward_time[b, a],
    forward_time[a, b], forward_time[b, a]
  )
  valid_time <- travel_time[is.finite(travel_time) & travel_time > 0]

  tibble(
    Sample_1 = a,
    Sample_2 = b,
    back_ab = connectivity[1],
    back_ba = connectivity[2],
    fwd_ab = connectivity[3],
    fwd_ba = connectivity[4],
    bmin_ab = travel_time[1],
    bmin_ba = travel_time[2],
    fmin_ab = travel_time[3],
    fmin_ba = travel_time[4],
    bi_connectivity_mean = mean(connectivity),
    bi_connectivity_max = max(connectivity),
    bi_connected = any(connectivity > 0),
    bi_min_travel_days = ifelse(length(valid_time) > 0, min(valid_time), NA_real_),
    directed_link_count = sum(connectivity > 0)
  )
}) |>
  bind_rows() |>
  mutate(pair_key = paste(pmin(Sample_1, Sample_2), pmax(Sample_1, Sample_2), sep = "__"))

assembly <- read_csv("output/raw_tables/assembly_module_pairwise.csv", show_col_types = FALSE) |>
  mutate(pair_key = paste(pmin(Sample_1, Sample_2), pmax(Sample_1, Sample_2), sep = "__"))

tina <- read_csv("output/raw_tables/tina_pina_pairwise.csv", show_col_types = FALSE) |>
  mutate(pair_key = paste(pmin(sample_1, sample_2), pmax(sample_1, sample_2), sep = "__")) |>
  select(pair_key, tina_dist, pina_dist, bray)

pairwise_data <- assembly |>
  left_join(ocean_pairs |> select(-Sample_1, -Sample_2), by = "pair_key") |>
  left_join(tina, by = "pair_key") |>
  mutate(
    bray_similarity = 1 - bray,
    tina_similarity = 1 - tina_dist
  )

write_csv(
  pairwise_data,
  "output/raw_tables/oceanparcels_stationdates_bidirectional_pairwise_365d50km7.csv"
)

pair_matrix <- function(value, diagonal = 0) {
  m <- matrix(diagonal, length(stations), length(stations), dimnames = list(stations, stations))
  for (i in seq_len(nrow(pairwise_data))) {
    a <- pairwise_data$Sample_1[i]
    b <- pairwise_data$Sample_2[i]
    m[a, b] <- pairwise_data[[value]][i]
    m[b, a] <- pairwise_data[[value]][i]
  }
  m
}

qap_cor <- function(y, x, mask, permutations = 9999, seed = 123456789) {
  idx <- which(lower.tri(x) & mask, arr.ind = TRUE)
  observed <- cor(y[idx], x[idx], method = "spearman")
  set.seed(seed)
  null <- replicate(permutations, {
    p <- sample(seq_len(nrow(y)))
    yp <- y[p, p]
    cor(yp[idx], x[idx], method = "spearman")
  })
  tibble(
    observed_r = observed,
    qap_p_two_sided = (sum(abs(null) >= abs(observed)) + 1) / (permutations + 1)
  )
}

module_matrix <- pair_matrix("module_similarity", diagonal = 1)
bray_matrix <- pair_matrix("bray_similarity", diagonal = 1)
travel_matrix <- pair_matrix("bi_min_travel_days", diagonal = NA_real_)
connectivity_matrix <- pair_matrix("bi_connectivity_mean", diagonal = 0)
connected_mask <- is.finite(travel_matrix)
all_pairs_mask <- matrix(TRUE, length(stations), length(stations))

qap_statistics <- bind_rows(
  qap_cor(module_matrix, travel_matrix, connected_mask) |>
    mutate(test = "module_similarity_vs_min_travel_days", .before = 1),
  qap_cor(module_matrix, connectivity_matrix, all_pairs_mask) |>
    mutate(test = "module_similarity_vs_connectivity_mean", .before = 1),
  qap_cor(bray_matrix, travel_matrix, connected_mask) |>
    mutate(test = "bray_similarity_vs_min_travel_days", .before = 1)
)
write_csv(qap_statistics, "output/raw_tables/oceanparcels_stationdates_qap_statistics.csv")

idx <- which(lower.tri(module_matrix), arr.ind = TRUE)
connected <- connectivity_matrix[idx] > 0
observed_difference <- median(module_matrix[idx][connected]) - median(module_matrix[idx][!connected])

set.seed(123456789)
null_difference <- replicate(9999, {
  p <- sample(seq_len(nrow(module_matrix)))
  permuted_module_matrix <- module_matrix[p, p]
  median(permuted_module_matrix[idx][connected]) - median(permuted_module_matrix[idx][!connected])
})

write_csv(
  tibble(
    test = "module_similarity_connected_minus_unconnected",
    median_difference = observed_difference,
    qap_p_two_sided = (sum(abs(null_difference) >= abs(observed_difference)) + 1) / 10000,
    permutations = 9999
  ),
  "output/raw_tables/oceanparcels_stationdates_connected_unconnected_qap.csv"
)

station_zones <- bind_rows(
  assembly |> select(station = Sample_1, oceanic_sector = Zone_1),
  assembly |> select(station = Sample_2, oceanic_sector = Zone_2)
) |>
  distinct()

station_connectivity <- bind_rows(
  pairwise_data |> select(station = Sample_1, bi_connected, bi_connectivity_mean),
  pairwise_data |> select(station = Sample_2, bi_connected, bi_connectivity_mean)
) |>
  group_by(station) |>
  summarise(
    link_count = sum(bi_connected),
    connectivity_sum = sum(bi_connectivity_mean),
    .groups = "drop"
  ) |>
  left_join(station_zones, by = "station") |>
  select(station, oceanic_sector, link_count, connectivity_sum)

write_csv(
  station_connectivity,
  "output/raw_tables/oceanparcels_stationdates_station_connectivity.csv"
)

station_connectivity |>
  group_by(oceanic_sector) |>
  summarise(
    n = n(),
    connectivity_mean = mean(connectivity_sum),
    connectivity_sd = sd(connectivity_sum),
    links_mean = mean(link_count),
    links_sd = sd(link_count),
    .groups = "drop"
  ) |>
  write_csv("output/raw_tables/oceanparcels_stationdates_zone_connectivity_summary.csv")
