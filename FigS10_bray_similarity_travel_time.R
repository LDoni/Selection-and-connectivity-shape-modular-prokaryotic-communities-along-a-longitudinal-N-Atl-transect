## FigS10  Bray similarity vs travel time

source("functions/01_atlantic_transect_and_community.R")
source("functions/04_oceanographic_connectivity.R")


 
pairs <- read_ocean_pairs() |>
  select(Sample_1, Sample_2, bi_min_travel_days, Mechanism) |>
  filter(!is.na(bi_min_travel_days)) |>
  mutate(pair_id = paste(pmin(Sample_1, Sample_2), pmax(Sample_1, Sample_2), sep = "__"))

bray_pairs <- read_csv("output/raw_tables/tina_pina_pairwise.csv", show_col_types = FALSE) |>
  transmute(pair_id = paste(pmin(sample_1, sample_2), pmax(sample_1, sample_2), sep = "__"), bray_similarity = 1 - bray)

pairs <- pairs |>
  left_join(bray_pairs, by = "pair_id")

ggplot(pairs, aes(bi_min_travel_days, bray_similarity)) +
  geom_point(color = "grey45", size = 2.5, alpha = 0.85) +
  geom_smooth(method = "lm", color = "black", fill = "grey75") +
  labs(x = "Combined forward-backward minimum travel time (days)", y = "Bray-Curtis similarity") +
  theme_bw()

cor.test(pairs$bi_min_travel_days, pairs$bray_similarity, method = "spearman", exact = FALSE)

#yes