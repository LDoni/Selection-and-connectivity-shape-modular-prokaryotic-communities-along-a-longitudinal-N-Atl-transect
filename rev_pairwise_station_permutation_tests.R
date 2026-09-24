#script for the revision of the pairwise station permutation with also fig 5B!

library(readr)
library(dplyr)
library(ggplot2)

# Input created in 04_prepare_ocean_connectivity_pairs.R
pairs <- read_csv(
  "output/raw_tables/oceanparcels_stationdates_bidirectional_pairwise_365d50km7.csv",
  show_col_types = FALSE
)

set.seed(123456789)
stations <- unique(c(pairs$Sample_1, pairs$Sample_2))
n <- length(stations)


#facciamo la func per la matrice
make_matrix <- function(x) {
  m <- matrix(NA_real_, n, n, dimnames = list(stations, stations))
  for (i in 1:nrow(pairs)) {
    a <- pairs$Sample_1[i]
    b <- pairs$Sample_2[i]
    m[a, b] <- x[i]
    m[b, a] <- x[i]
  }
  diag(m) <- NA
  m
}


module <- make_matrix(pairs$module_similarity)
bray <- make_matrix(pairs$bray_similarity)
connectivity <- make_matrix(pairs$bi_connectivity_mean)
travel <- make_matrix(pairs$bi_min_travel_days)
sst <- make_matrix(pairs$temp_diff)
longitude <- make_matrix(pairs$lon_diff)
connected <- make_matrix(as.numeric(pairs$bi_connected))



all_pairs <- upper.tri(module)
travel_pairs <- all_pairs & is.finite(travel)
connected_pairs <- all_pairs & connected == 1
unconnected_pairs <- all_pairs & connected == 0



# Spearman  with station-label permutations (9999 perm)
perm_cor <- function(y, x, mask) {
  observed <- cor(y[mask], x[mask], method = "spearman", use = "complete.obs")
  null <- numeric(9999)

  for (i in 1:9999) {
    p <- sample(1:n)
    y_perm <- y[p, p]
    null[i] <- cor(y_perm[mask], x[mask],
                   method = "spearman", use = "complete.obs")
  }

  p_value <- (sum(abs(null) >= abs(observed)) + 1) / (9999 + 1)
  c(statistic = observed, p = p_value)
}


# Difference medians with station-label permutations )9999 perm)
perm_median <- function(y, group1, group2) {
  observed <- median(y[group1], na.rm = TRUE) -
              median(y[group2], na.rm = TRUE)
  null <- numeric(9999)

  for (i in 1:9999) {
    p <- sample(1:n)
    y_perm <- y[p, p]
    null[i] <- median(y_perm[group1], na.rm = TRUE) -
               median(y_perm[group2], na.rm = TRUE)
  }

  p_value <- (sum(abs(null) >= abs(observed)) + 1) / (9999 + 1)
  c(statistic = observed, p = p_value)
}

## new test!!!!
perm_cor(module, connectivity, all_pairs
perm_cor(module, travel, travel_pairs)
perm_cor(bray, travel, travel_pairs)

perm_median(module, connected_pairs, unconnected_pairs)
perm_median(sst, connected_pairs, unconnected_pairs)
perm_median(longitude, connected_pairs, unconnected_pairs)


# Fig. 5B revised
fit_data <- pairs |> filter(is.finite(bi_min_travel_days))
grid <- seq(min(fit_data$bi_min_travel_days),
            max(fit_data$bi_min_travel_days), length.out = 100)

fit_all <- lm(module_similarity ~ bi_min_travel_days, data = fit_data)
observed_fit <- predict(fit_all,
                        newdata = data.frame(bi_min_travel_days = grid))

jack_fits <- matrix(NA_real_, nrow = n, ncol = length(grid))

for (i in 1:n) {
  station_out <- stations[i]
  keep <- fit_data$Sample_1 != station_out &
          fit_data$Sample_2 != station_out

  fit_i <- lm(module_similarity ~ bi_min_travel_days,
              data = fit_data[keep, ])

  jack_fits[i, ] <- predict(
    fit_i,
    newdata = data.frame(bi_min_travel_days = grid)
  )
}

jack_mean <- colMeans(jack_fits)
jack_se <- sqrt((n - 1) / n *
                colSums((jack_fits - jack_mean)^2))
t_value <- qt(0.975, df = n - 1)

plot_data <- data.frame(
  travel_time = grid,
  fit = observed_fit,
  lower = observed_fit - t_value * jack_se,
  upper = observed_fit + t_value * jack_se
)



         #plot with ggplot
  ggplot(fit_data, aes(bi_min_travel_days, module_similarity)) +
  geom_point() +
  geom_ribbon(data = plot_data,
              aes(x = travel_time, ymin = lower, ymax = upper),
              inherit.aes = FALSE, alpha = 0.25) +
  geom_line(data = plot_data,
            aes(x = travel_time, y = fit),
            inherit.aes = FALSE) +
  labs(x = "Minimum travel time (days)",
       y = "Module similarity") +
  theme_classic()








         
         

