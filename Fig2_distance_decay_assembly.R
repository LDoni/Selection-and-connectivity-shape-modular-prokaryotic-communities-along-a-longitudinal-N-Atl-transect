## Figure 2 - distance decay and assembly processes

source("functions/01_atlantic_transect_and_community.R")
source("functions/02_assembly_processes.R")

ps <- read_atlantic_ps()
meta <- sample_metadata(ps)

bray <- distance(ps, method = "bray")
geo <- as.dist(geo_matrix(meta))
temp <- dist(meta$sst_degC)

dd <- data.frame(
  geo_km = as.vector(geo),
  temp_diff = as.vector(temp),
  bray_similarity = 1 - as.vector(bray)
)


### Fig 2A (distance decay bray sim)

dd_long <- dd |>
  transmute(
    bray_similarity,
    `Geographic distance` = rescale(geo_km),
    `Temperature difference` = rescale(temp_diff)
  ) |>
  pivot_longer(-bray_similarity, names_to = "Distance", values_to = "Scaled_distance")

ggplot(dd_long, aes(Scaled_distance, bray_similarity, color = Distance, linetype = Distance)) +
  geom_point(alpha = 0.55, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("Geographic distance" = "#2166ac", "Temperature difference" = "#b2182b")) +
  scale_linetype_manual(values = c("Geographic distance" = "solid", "Temperature difference" = "dashed")) +
  labs(x = "Scaled distance", y = "Bray-Curtis similarity") +
  theme_bw()

bray_geo_lm <- lm(bray_similarity ~ geo_km, data = dd)
bray_temp_lm <- lm(bray_similarity ~ temp_diff, data = dd)
set.seed(123456789)
bray_geo_mantel <- mantel(bray, geo, permutations = 9999)
set.seed(123456789)
bray_geo_partial_mantel_temperature <- mantel.partial(bray, geo, temp, permutations = 9999)
set.seed(123456789)
bray_temp_mantel <- mantel(bray, temp, permutations = 9999)

summary(bray_geo_lm)
summary(bray_temp_lm)
bray_geo_mantel
bray_geo_partial_mantel_temperature
bray_temp_mantel


### Fig 2B (ecological mechanisms)

assembly <- read_assembly_pairs()

set.seed(123456789)
boot_mech <- replicate(
  5000,
  {
    x <- sample(assembly$Mechanism, replace = TRUE)
    prop.table(table(factor(x, levels = unique(assembly$Mechanism))))
  }
)

global <- assembly |>
  count(Mechanism) |>
  mutate(prop = n / sum(n)) |>
  left_join(
    data.frame(
      Mechanism = rownames(boot_mech),
      low = apply(boot_mech, 1, quantile, 0.025, na.rm = TRUE),
      high = apply(boot_mech, 1, quantile, 0.975, na.rm = TRUE)
    ),
    by = "Mechanism"
  )

ggplot(global, aes(reorder(Mechanism, prop), prop, fill = Mechanism)) +
  geom_col(color = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.18) +
  coord_flip() +
  scale_y_continuous(labels = percent_format()) +
  labs(x = NULL, y = "Pairwise comparisons") +
  theme_bw() +
  theme(legend.position = "none")


### Fig 2C(ecological mechanisms distance)

assembly$distance_bin <- cut(
  assembly$geographic_km,
  c(0, 250, 500, 1000, 2000, 3000, 5000, Inf),
  include.lowest = TRUE,
  right = FALSE
)

by_distance <- assembly |>
  count(distance_bin, Mechanism) |>
  group_by(distance_bin) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

ggplot(by_distance, aes(distance_bin, prop, fill = Mechanism)) +
  geom_col(color = "black", linewidth = 0.15) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Geographic distance (km)", y = "Pairwise comparisons") +
  theme_bw()


### Fig 2D (ecological mechanisms Temp )

assembly$temp_bin <- cut(assembly$temp_diff, c(0, 1, 2, 4, 6, 8, Inf), include.lowest = TRUE, right = FALSE)

by_temp <- assembly |>
  count(temp_bin, Mechanism) |>
  group_by(temp_bin) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

ggplot(by_temp, aes(temp_bin, prop, fill = Mechanism)) +
  geom_col(color = "black", linewidth = 0.15) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Temperature difference class (deg C)", y = "Pairwise comparisons") +
  theme_bw()
