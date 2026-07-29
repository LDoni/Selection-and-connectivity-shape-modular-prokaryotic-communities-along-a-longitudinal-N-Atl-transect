## Figure S9 - variation partitioning

source("functions/01_atlantic_transect_and_community.R")
source("functions/03_network_and_modules.R")

### Fig S9

vp <- read_csv("output/raw_tables/variation_partitioning_summary.csv", show_col_types = FALSE) |>
  mutate(
    Fraction = recode(Fraction,
      `Space|Environment` = "Space",
      `Environment|Space` = "Environment",
      `Environment+Space` = "Environment+Space"),
    Fraction = factor(Fraction, levels = c("Environment+Space", "Environment", "Space")),
    Response = factor(Response, levels = c("Community", "Main modules"))
  )

ggplot(vp, aes(Adj_R2, Fraction, fill = Fraction)) +
  geom_col(width = 0.68, color = "black", linewidth = 0.25) +
  facet_wrap(~Response, nrow = 1) +
  scale_fill_manual(values = c("Space" = "brown", "Environment" = "gray",
                               "Environment+Space" = "green"), guide = "none") +
  labs(x = expression(Adjusted~R^2), y = NULL) +
  theme_bw()

