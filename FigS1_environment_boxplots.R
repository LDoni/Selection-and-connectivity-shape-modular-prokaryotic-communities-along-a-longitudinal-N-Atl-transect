## Figure S1 - environmental boxplots

source("functions/01_atlantic_transect_and_community.R")

ps <- read_atlantic_ps()
meta <- sample_metadata(ps)


### Fig S1 env variables

env <- meta |>
  select(sample_id, oceanic_sector, sst_degC, salinity_psu, chlorophyll_mg_m3, oxygen_mmol_m3, nitrate_mmol_m3, phosphate_mmol_m3, silicate_mmol_m3) |>
  pivot_longer(-c(sample_id, oceanic_sector), names_to = "Variable", values_to = "Value") |>
  mutate(Variable = recode(Variable,
    sst_degC = "Temperature (°C)",
    salinity_psu = "Salinity (PSU)",
    chlorophyll_mg_m3 = "Chlorophyll a (mg m-3)",
    oxygen_mmol_m3 = "Dissolved oxygen (mmol m-3)",
    nitrate_mmol_m3 = "Nitrate (mmol m-3)",
    phosphate_mmol_m3 = "Phosphate (mmol m-3)",
    silicate_mmol_m3 = "Silicate (mmol m-3)"
  ))

ggplot(env, aes(oceanic_sector, Value, fill = oceanic_sector, color = oceanic_sector)) +
  geom_boxplot(width = 0.62, alpha = 0.35, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.10), size = 2) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = zone_cols, drop = FALSE) +
  scale_color_manual(values = zone_cols, drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = "none")

