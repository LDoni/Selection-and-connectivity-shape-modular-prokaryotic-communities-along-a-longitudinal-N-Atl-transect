## Figure S7 - module taxonomy

library(readr)
library(dplyr)
library(ggplot2)
library(scales)


### Fig S7 taxonomy main modules

tax_family <- read_csv(
  "output/raw_tables/main_module_taxonomy_composition_all_ranks.csv",
  show_col_types = FALSE
) |>
  filter(Rank == "Family", Module %in% 1:4)

top_families <- tax_family |>
  group_by(Taxon) |>
  summarise(total = sum(Percent_abundance_within_module), .groups = "drop") |>
  slice_max(total, n = 10) |>
  pull(Taxon)

tax_family <- tax_family |>
  mutate(Family = ifelse(Taxon %in% top_families, Taxon, "Other")) |>
  group_by(Module, Family) |>
  summarise(percent = sum(Percent_abundance_within_module), .groups = "drop")

ggplot(tax_family, aes(paste0("Module ", Module), percent, fill = Family)) +
  geom_col(color = "white", linewidth = 0.1) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Family composition (%)") +
  theme_bw()

