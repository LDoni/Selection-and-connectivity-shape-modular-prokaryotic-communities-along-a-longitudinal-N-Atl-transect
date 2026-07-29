## Figure 1 - community structure

source("functions/01_atlantic_transect_and_community.R")

ps <- read_atlantic_ps()
 
meta <- sample_metadata(ps)

map_df <- meta |>
  mutate(Longitude = longitude, Latitude = latitude)

#FIG1A (transect map w ann temp)

# 2022 annual mean sea surface temperatures; source: https://oceandata.sci.gsfc.nasa.gov/file_search/
# to download the file: system('wget --auth-no-challenge --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies "https://oceandata.sci.gsfc.nasa.gov/ob/getfile/AQUA_MODIS.20220101_20221231.L3m.YR.SST.sst.4km.nc"')

library(RColorBrewer)

sst_file <- Sys.getenv(
  "ATLANTIC_SST_FILE",
  unset = "AQUA_MODIS.20220101_20221231.L3m.YR.SST.sst.4km.nc"
)

world <- map_data("world")

p_map <- ggplot() +
  geom_polygon(
    data = world,
    aes(long, lat, group = group),
    fill = "#D9D9D9",
    color = "gray70",
    linewidth = 0.2
  )


  nc <- ncdf4::nc_open(sst_file)
  lon <- ncdf4::ncvar_get(nc, "lon")
  lat <- ncdf4::ncvar_get(nc, "lat")
  lon_i <- which(lon >= -90 & lon <= 10)
  lat_i <- which(lat >= 25 & lat <= 45)
  sst <- ncdf4::ncvar_get(
    nc, "sst",
    start = c(min(lon_i), min(lat_i)),
    count = c(length(lon_i), length(lat_i))
  )
  ncdf4::nc_close(nc)

  SST_df <- expand.grid(x = lon[lon_i], y = lat[lat_i])
  SST_df$sst <- as.vector(sst)

  p_map <- ggplot() +
    geom_raster(data = SST_df, aes(x = x, y = y, fill = sst)) +
    scale_fill_gradientn(
      colors = rev(brewer.pal(11, "RdBu")),
      name = "Annual SST (°C)",
      limits = c(0, 35),
      na.value = "transparent"
    ) +
    geom_polygon(
      data = world,
      aes(long, lat, group = group),
      fill = "#D9D9D9",
      color = "gray70",
      linewidth = 0.2
    )


p_map +
  geom_point(data = map_df, aes(Longitude, Latitude, color = oceanic_sector), size = 3) +
  geom_text(data = map_df, aes(Longitude, Latitude, label = day_of_navigation), size = 3, nudge_y = 1, check_overlap = TRUE) +
  scale_color_manual(values = zone_cols, drop = FALSE) +
  coord_quickmap(xlim = c(-90, 10), ylim = c(25, 45), expand = FALSE) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)", color = "Oceanic area") +
  theme_minimal()


### Fig 1B/C (richness)


pa <- transform_sample_counts(ps, function(x) ifelse(x > 0, 1, 0))
rich <- estimate_richness(pa, measures = "Observed") |>
  rownames_to_column("sample_id") |>
  left_join(meta |> select(sample_id, oceanic_sector, longitude), by = "sample_id")
# U plot
ggplot(rich, aes(longitude, Observed, color = oceanic_sector)) +
  geom_point(size = 3) +
  geom_smooth(color = "black", se = TRUE) +
  scale_color_manual(values = zone_cols, drop = FALSE) +
  labs(x = "Longitude (°E)", y = "Observed ASVs", color = "Oceanic area") +
  theme_bw()

#boxplots
ggplot(rich, aes(oceanic_sector, Observed, fill = oceanic_sector)) +
  geom_boxplot(alpha = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2.4) +
  scale_fill_manual(values = zone_cols, drop = FALSE) +
  labs(x = NULL, y = "Observed ASVs") +
  theme_bw() +
  theme(legend.position = "none")



# Fig 1D (PCoA)

library(ggside)
library(ggforce)
library(ggrepel)

bray <- distance(ps, method = "bray")
ord_bray <- ordinate(ps, method = "PCoA", distance = bray)

plot_ordination(ps, ord_bray,  color="oceanic_sector",axes =c(1,2))+
  scale_color_manual(values = c("#4280fc", "#ffb452","#f7170a"))+ 
  geom_point(size=3)+ 
  geom_text_repel(aes(label=day_of_navigation),max.overlaps = Inf, show.legend = FALSE)+
  geom_mark_ellipse(aes(color = oceanic_sector), show.legend = FALSE)+ 
  theme_void()+theme_bw()

set.seed(123456789)
adonis2(bray ~ oceanic_sector, data = meta)
anova(betadisper(bray, meta$oceanic_sector))


# Fig 1E (comm composition barplot)


ps_comp <- read_atlantic_raw_ps()
ps_comp <- prune_taxa(taxa_names(ps), ps_comp)
ps_comp <- tss_normalize_ps(ps_comp)

otu_comp <- as(otu_table(ps_comp), "matrix")
if (!taxa_are_rows(ps_comp)) {
  otu_comp <- t(otu_comp)
}

tax_comp <- as.data.frame(tax_table(ps_comp), stringsAsFactors = FALSE)
class_comp <- tax_comp$Class
class_comp[is.na(class_comp) | class_comp == ""] <- "Unclassified"

class_mat <- rowsum(otu_comp, class_comp)
class_avg <- sort(rowMeans(class_mat), decreasing = TRUE)
keep_classes <- names(class_avg)[names(class_avg) != "Unclassified"][1:10]
class_levels <- c(keep_classes, "Others")

class_mat <- rbind(
  class_mat[keep_classes, , drop = FALSE],
  Others = colSums(class_mat[setdiff(rownames(class_mat), keep_classes), , drop = FALSE])
)

class_meta <- sample_metadata(ps_comp) |>
  transmute(
    sample_id,
    station = sample_id,
    oceanic_sector,
    longitude
  )

class_df <- as.data.frame(class_mat) |>
  rownames_to_column("Class") |>
  pivot_longer(-Class, names_to = "sample_id", values_to = "Relative_abundance") |>
  left_join(class_meta, by = "sample_id") |>
  mutate(
    station = factor(station, levels = unique(station[order(longitude)])),
    oceanic_sector = factor(oceanic_sector, levels = names(zone_cols)),
    Class = factor(Class, levels = class_levels)
  )

ggplot(class_df, aes(station, Relative_abundance, fill = Class)) +
  geom_col(width = 0.9) +
  facet_grid(~ oceanic_sector, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = NULL, y = "Relative abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
