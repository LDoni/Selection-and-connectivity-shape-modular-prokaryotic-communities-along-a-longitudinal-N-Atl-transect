## Figure 4 - module environment

source("functions/01_atlantic_transect_and_community.R")
source("functions/03_network_and_modules.R")

ps <- read_atlantic_ps()
meta <- sample_metadata(ps)

mod_mat <- module_abundance(ps)
main_mod <- top_modules(mod_mat, 4)

corr <- read_csv("output/raw_tables/module_environment_correlations.csv", show_col_types = FALSE)
labs <- c(Longitude.x = "Longitude", so = "Salinity", chl = "Chlorophyll",
          o2 = "Oxygen", no3 = "Nitrate", po4 = "Phosphate",
          si = "Silicate", thetao = "Temperature")


### Fig 4A (main modules correlogram)

corr_plot <- corr |>
  filter(Module %in% main_mod, Variable %in% names(labs)) |>
  mutate(
    Module = factor(paste0("Module ", Module), levels = paste0("Module ", main_mod)),
    Variable = factor(labs[Variable], levels = unname(labs)),
    stars = case_when(fdr < 0.001 ~ "***", fdr < 0.01 ~ "**", fdr < 0.05 ~ "*", TRUE ~ "")
  )

ggplot(corr_plot, aes(Module, Variable, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = stars), size = 4) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, fill = "Spearman r") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### Fig 4B termal optimum and breadth

niche <- read_csv("output/raw_tables/main_module_thermal_niche.csv", show_col_types = FALSE) |>
  mutate(Module = factor(Module, levels = c(1, 4, 2, 3)))

ggplot(niche, aes(thermal_optimum, Module)) +
  geom_segment(aes(x = q10, xend = q90, yend = Module), color = "grey70", linewidth = 1.2) +
  geom_segment(aes(x = thermal_optimum - thermal_breadth_sd,
                   xend = thermal_optimum + thermal_breadth_sd,
                   yend = Module, color = thermal_optimum), linewidth = 4) +
  geom_point(aes(size = mean_relative_abundance, fill = thermal_optimum), shape = 21, color = "black") +
  geom_text(aes(label = Module), color = "black", size = 3.6) +
  scale_color_gradient2(low = "#2b83ba", mid = "#abdda4", high = "#f46d43", midpoint = 23) +
  scale_fill_gradient2(low = "#2b83ba", mid = "#abdda4", high = "#f46d43", midpoint = 23,
                       name = "Thermal optimum") +
  scale_size_continuous(range = c(5, 13), breaks = c(0.10, 0.25, 0.45),
                        labels = percent, name = "Mean relative abundance") +
  scale_x_continuous(breaks = c(20, 22, 24, 26, 28), expand = expansion(mult = c(0.02, 0.16))) +
  scale_y_discrete(limits = c("1", "4", "2", "3")) +
  labs(x = "Thermal optimum and breadth (°C)", y = "Dominant module") +
  guides(color = "none") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

