## Figure S2 - weighted UniFrac PCoA

source("functions/01_atlantic_transect_and_community.R")

ps <- read_atlantic_ps()


### Fig S2 pcoa wunifrac
library(ggrepel)
wuni <- distance(ps, method = "wunifrac")
ord <- ordinate(ps, method = "PCoA", distance = wuni)



plot_ordination(ps, ord, color = "oceanic_sector") +
  geom_point(size = 3) +
  scale_color_manual(values = zone_cols, drop = FALSE) +
  labs(
    x = paste0("PCoA1 (", round(100 * ord$values$Relative_eig[1], 1), "%)"),
    y = paste0("PCoA2 (", round(100 * ord$values$Relative_eig[2], 1), "%)"),
    color = "Oceanic area"
  ) +
  theme_bw()+
  geom_text_repel(
    aes(label = sample_id ),
    size = 3,
    max.overlaps = Inf)

meta <- sample_metadata(ps)
set.seed(123456789)
permanova_wunifrac_zone <- adonis2(wuni ~ oceanic_sector, data = meta)
betadisper_wunifrac_zone <- anova(betadisper(wuni, meta$oceanic_sector))

permanova_wunifrac_zone
betadisper_wunifrac_zone
