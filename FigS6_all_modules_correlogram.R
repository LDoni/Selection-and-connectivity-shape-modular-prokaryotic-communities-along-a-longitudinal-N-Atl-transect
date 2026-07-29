## Figure S6 - all modules correlogram

source("functions/01_atlantic_transect_and_community.R")
source("functions/03_network_and_modules.R")

ps <- read_atlantic_ps()
meta <- sample_metadata(ps)
mod_mat <- module_abundance(ps)




corr <- module_correlations(mod_mat, meta)
labs <- c(sst_degC = "SST", salinity_psu = "Salinity", chlorophyll_mg_m3 = "Chlorophyll", oxygen_mmol_m3 = "Oxygen", nitrate_mmol_m3 = "Nitrate", phosphate_mmol_m3 = "Phosphate", silicate_mmol_m3 = "Silicate", longitude = "Longitude", latitude = "Latitude")

corr_plot <- corr |>
  mutate(
    Module = factor(paste0("Module ", Module), levels = paste0("Module ", rownames(mod_mat))),
    Variable = factor(labs[Variable], levels = labs),
    stars = case_when(FDR < 0.001 ~ "***", FDR < 0.01 ~ "**", FDR < 0.05 ~ "*", TRUE ~ "")
  )

ggplot(corr_plot, aes(Module, Variable, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = stars), size = 3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, fill = "Spearman r") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

