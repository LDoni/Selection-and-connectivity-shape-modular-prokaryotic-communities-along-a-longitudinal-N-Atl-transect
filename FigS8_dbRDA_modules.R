

source("functions/01_atlantic_transect_and_community.R")
source("functions/03_network_and_modules.R")
library(ggrepel)

ps <- read_atlantic_ps()
meta <- sample_metadata(ps)
mod_mat <- module_abundance(ps)
main_mod <- top_modules(mod_mat, 4)


### Fig S8dbRDA modules

comm <- t(mod_mat[main_mod, , drop = FALSE])
comm_hel <- decostand(comm, "hellinger")
env <- meta[match(rownames(comm_hel), meta$sample_id), c("sst_degC", "oxygen_mmol_m3", "nitrate_mmol_m3", "phosphate_mmol_m3", "silicate_mmol_m3", "salinity_psu", "chlorophyll_mg_m3", "longitude", "latitude")]
env <- data.frame(lapply(env, as.numeric), row.names = rownames(comm_hel))

fit <- capscale(comm_hel ~ sst_degC + oxygen_mmol_m3 + nitrate_mmol_m3 + phosphate_mmol_m3 + silicate_mmol_m3 + salinity_psu + chlorophyll_mg_m3 + longitude + latitude, data = env, distance = "euclidean")
scr <- scores(fit, display = "sites", scaling = 2)
scr <- data.frame(scr, sample_id = rownames(scr)) |>
  left_join(meta |> select(sample_id, oceanic_sector), by = "sample_id")
mods <- data.frame(scores(fit, display = "species", scaling = 2)) |>
  rownames_to_column("Module") |>
  mutate(Module = sub("^X", "", Module))
env_scores <- data.frame(scores(fit, display = "bp", scaling = 2)) |>
  rownames_to_column("Variable") |>
  mutate(Variable = recode(Variable,
    sst_degC = "Temperature", oxygen_mmol_m3 = "Oxygen", nitrate_mmol_m3 = "Nitrate",
    phosphate_mmol_m3 = "Phosphate", silicate_mmol_m3 = "Silicate",
    salinity_psu = "Salinity", chlorophyll_mg_m3 = "Chlorophyll",
    longitude = "Longitude", latitude = "Latitude"))

ggplot(scr, aes(CAP1, CAP2, color = oceanic_sector)) +
  geom_point(size = 3) +
  stat_ellipse(linewidth = 0.7, show.legend = FALSE) +
  geom_segment(data = env_scores, aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               inherit.aes = FALSE, arrow = arrow(length = grid::unit(0.16, "cm")), color = "grey25") +
  geom_text_repel(data = env_scores, aes(CAP1, CAP2, label = Variable),
                  inherit.aes = FALSE, color = "grey25", size = 3) +
  geom_text_repel(data = mods, aes(CAP1, CAP2, label = Module),
                  inherit.aes = FALSE, color = "black", fontface = "bold", size = 4) +
  scale_color_manual(values = zone_cols, drop = FALSE) +
  labs(x = "dbRDA1", y = "dbRDA2", color = "Oceanic area") +
  theme_bw()

set.seed(123456789)
anova.cca(fit, by = "terms", permutations = 999)
set.seed(123456789)
anova.cca(fit, by = "axis", permutations = 999)
