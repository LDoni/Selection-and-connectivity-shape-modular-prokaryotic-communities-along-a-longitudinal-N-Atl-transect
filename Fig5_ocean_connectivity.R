## Figure 5 - oceanographic connectivity

source("functions/04_oceanographic_connectivity.R")
library(ggplot2)


### Fig 5A


## Figure 5A - oceanographic connectivity map

source("functions/04_oceanographic_connectivity.R")
library(ggplot2)
library(dplyr)
library(readr)
library(maps)

pairs <- read_ocean_pairs()
stations <- read_csv("data/stations.csv", show_col_types = FALSE)

station_data <- bind_rows(
  pairs |> select(station = Sample_1, connectivity = bi_connectivity_mean),
  pairs |> select(station = Sample_2, connectivity = bi_connectivity_mean)
) |>
  group_by(station) |>
  summarise(connectivity = sum(connectivity, na.rm = TRUE), .groups = "drop") |>
  left_join(stations, by = "station")

links <- pairs |>
  filter(bi_connected, !is.na(bi_min_travel_days)) |>
  select(Sample_1, Sample_2, bi_connectivity_mean, bi_min_travel_days) |>
  left_join(stations |> select(Sample_1 = station, lat_1 = lat, lon_1 = lon), by = "Sample_1") |>
  left_join(stations |> select(Sample_2 = station, lat_2 = lat, lon_2 = lon), by = "Sample_2")

world <- map_data("world") |>
  filter(
    long > -90,
    long < 10,
    lat > 15,
    lat < 50
  )

ggplot() +
  geom_polygon( data = world,  aes( x = long,
                                    y = lat,  group = group),
                fill = "grey88",color = "grey45", linewidth = 0.25
  ) +
  geom_segment(data = links,aes(x = lon_1,
                                y = lat_1,  xend = lon_2,  yend = lat_2,color = bi_min_travel_days,
                                linewidth = bi_connectivity_mean    ),
               alpha = 0.9,lineend = "round"  ) +
  geom_point(  data = station_data,
               aes(x = lon,y = lat,size = connectivity),
               shape = 21,fill = "white",color = "black",stroke = 0.35) +
  scale_color_gradient(
    low = "#d73027",
    high = "#fff7ec",
    name = "Minimum travel time (days)"
  ) +
  scale_linewidth_continuous(
    range = c(0.35, 1.6),
    name = "Mean bidirectional connectivity"
  ) +
  scale_size_continuous(
    range = c(2.5, 7),
    name = "Station connectivity"
  ) +
  coord_quickmap( xlim = c(-83, -3),
                  ylim = c(23, 43), expand = FALSE) +
  labs(  x = "Longitude", y = "Latitude") +
  theme_bw() +  theme(panel.background = element_rect(fill = "lightblue", #background used in the real plot is too heavy
                                                      color = NA),panel.grid = element_blank())







### Connectivity by oceanic sector (not included)

# station_data <- read_csv(
#   "output/raw_tables/oceanparcels_stationdates_station_connectivity.csv",
#   show_col_types = FALSE
# ) |>
#   mutate(oceanic_sector = factor(oceanic_sector, levels = c("ANW", "ANC", "ANE")))
#
# ggplot(station_data, aes(oceanic_sector, connectivity_sum, fill = oceanic_sector)) +
#   geom_boxplot(width = 0.52, alpha = 0.68, outlier.shape = NA) +
#   geom_jitter(shape = 21, size = 3.2, width = 0.10) +
#   scale_fill_manual(values = c(ANW = "#4280fc", ANC = "#ffb452", ANE = "#f7170a"), guide = "none") +
#   labs(x = "Oceanic sector", y = "Station-level summed connectivity") +
#   theme_classic()


### Fig 5B

pairs <- read_ocean_pairs()

plot_df <- pairs |>
  filter(!is.na(bi_min_travel_days), !is.na(module_similarity))

ggplot(plot_df, aes(bi_min_travel_days, module_similarity)) +
  geom_point(color = "grey45", size = 2.5, alpha = 0.85) +
  geom_smooth(method = "lm", color = "black", fill = "grey75") +
  labs(x = "travel time (days)", y = "Module similarity") +
  theme_bw()

cor.test(plot_df$bi_min_travel_days, plot_df$module_similarity,
         method = "spearman", exact = FALSE)
summary(lm(module_similarity ~ bi_min_travel_days, data = plot_df))

