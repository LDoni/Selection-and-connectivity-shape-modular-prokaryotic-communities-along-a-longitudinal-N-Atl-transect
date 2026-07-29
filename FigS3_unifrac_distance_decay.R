## Figure S3 - weighted UniFrac distance decay

source("functions/01_atlantic_transect_and_community.R")

ps <- read_atlantic_ps()
meta <- sample_metadata(ps)


### Fig S3

uni <- distance(ps, method = "unifrac", weighted = TRUE)
geo <- as.dist(geo_matrix(meta))
temp <- dist(meta$sst_degC)

dd <- data.frame(
  geo_km = as.vector(geo),
  temp_diff = as.vector(temp),
  unifrac_similarity = 1 - as.vector(uni)
)

dd_long <- dd |>
  transmute(
    unifrac_similarity,
    `Geographic distance` = rescale(geo_km),
    `Temperature difference` = rescale(temp_diff)
  ) |>
  pivot_longer(-unifrac_similarity, names_to = "Distance", values_to = "Scaled_distance")

ggplot(dd_long, aes(Scaled_distance, unifrac_similarity, color = Distance, linetype = Distance)) +
  geom_point(alpha = 0.55, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("Geographic distance" = "#2166ac", "Temperature difference" = "#b2182b")) +
  scale_linetype_manual(values = c("Geographic distance" = "solid", "Temperature difference" = "dashed")) +
  labs(x = "Scaled distance", y = "Weighted UniFrac similarity") +
  theme_bw()

set.seed(123456789)
unifrac_geo_mantel <- mantel(uni, geo, permutations = 9999)
set.seed(123456789)
unifrac_geo_partial_mantel_temperature <- mantel.partial(uni, geo, temp, permutations = 9999)
set.seed(123456789)
unifrac_temp_mantel <- mantel(uni, temp, permutations = 9999)

unifrac_geo_mantel
unifrac_geo_partial_mantel_temperature
unifrac_temp_mantel
