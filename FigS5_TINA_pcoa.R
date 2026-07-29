## Figure S5 - TINA PCoA

source("functions/01_atlantic_transect_and_community.R")

pairwise <- read_csv("output/raw_tables/tina_pina_pairwise.csv", show_col_types = FALSE)
ps <- read_atlantic_ps()
meta <- sample_metadata(ps)


### Fig S5 TINA PCoA

samples <- meta$sample_id
tina_mat <- matrix(0, length(samples), length(samples), dimnames = list(samples, samples))
for (i in seq_len(nrow(pairwise))) {
  tina_mat[pairwise$sample_1[i], pairwise$sample_2[i]] <- pairwise$tina_dist[i]
  tina_mat[pairwise$sample_2[i], pairwise$sample_1[i]] <- pairwise$tina_dist[i]
}

ord <- cmdscale(as.dist(tina_mat), eig = TRUE, k = 2)
eig_pos <- ord$eig[ord$eig > 0]
var_exp <- round(100 * ord$eig[1:2] / sum(eig_pos), 1)

ord_df <- data.frame(sample_id = rownames(ord$points), PCoA1 = ord$points[, 1], PCoA2 = ord$points[, 2]) |>
  left_join(meta, by = "sample_id")

ggplot(ord_df, aes(PCoA1, PCoA2)) +
  geom_point(aes(fill = oceanic_sector), shape = 21, color = "black", size = 5.6) +
  ggrepel::geom_label_repel(aes(label = day_of_navigation), color = "black", size = 3) +
  scale_fill_manual(values = zone_cols, drop = FALSE) +
  coord_equal() +
  labs(x = paste0("PCoA1 (", var_exp[1], "%)"), y = paste0("PCoA2 (", var_exp[2], "%)"), fill = "Oceanic area") +
  theme_bw()

