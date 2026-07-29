## Figure 3 - network and modules

source("functions/01_atlantic_transect_and_community.R")
source("functions/03_network_and_modules.R")

ps <- read_atlantic_ps()
meta <- sample_metadata(ps)

g <- read_graph("output/raw_tables/SparCC_Network_pos.graphml", format = "graphml")
cl <- readRDS("output/raw_tables/cluster.rds")
cl <- cl[match(V(g)$name, cl$OTU_ID), ]
lay <- network_layout(g)
v_size <- network_node_size(g)
module_cols <- setNames(brewer.pal(12, "Set3"), as.character(1:12))


### Fig 3A --> network modules plot 

plot(g, layout = lay, vertex.size = v_size, vertex.label = NA,
     vertex.color = module_cols[as.character(cl$Cluster)],
     edge.color = adjustcolor("grey70", alpha.f = 0.25), main = "SparCC modules")


### Fig 3Bnetwork modules plot colored with Temp

plot(g, layout = lay, vertex.size = v_size, vertex.label = NA,
     vertex.color = network_continuous_colors(V(g)$sst_degC, "RdBu", reverse = TRUE),
     edge.color = adjustcolor("grey70", alpha.f = 0.25), main = "Temperature")
network_continuous_scale(V(g)$sst_degC, "RdBu", reverse = TRUE, label = "°C")


### Fig 3C--> network modules plot colored with longitude

plot(g, layout = lay, vertex.size = v_size, vertex.label = NA,
     vertex.color = network_continuous_colors(V(g)$longitude, "YlOrBr"),
     edge.color = adjustcolor("grey70", alpha.f = 0.25), main = "Longitude")
network_continuous_scale(V(g)$longitude, "YlOrBr", label = "°E")


## Fig 3D (modularity vs longitude)

net <- read.csv("output/raw_tables/Network_Properties.csv", check.names = FALSE) |>
  transmute(sample_id = Sample_ID, Modularity) |>
  left_join(meta |> select(sample_id, longitude, oceanic_sector), by = "sample_id")

ggplot(net, aes(longitude, 
                Modularity, 
                color = oceanic_sector)) +
  geom_point(size = 3) +
  geom_smooth(method = "gam", 
               color = "black") +
  scale_color_manual(values = zone_cols, drop = FALSE) +
  labs(x = "Longitude (°E)", y = "Local modularity") +
  theme_bw()

#not included in the paper
ggplot(net, aes(longitude, 
                Modularity, 
                color = oceanic_sector)) +
  geom_point(size = 3) +
  geom_smooth(method = "gam", 
              color = "black") +
  geom_smooth(method = "lm", 
              color = "red") +
  scale_color_manual(values = zone_cols, drop = FALSE) +
  labs(x = "Longitude (°E)", y = "Local modularity") +
  theme_bw()





local_modularity_longitude_gam <- mgcv::gam(Modularity ~ s(longitude, bs = "cs"), data = net, method = "REML")
summary(local_modularity_longitude_gam)


### Fig 3E Main modules (1-4) shift along long

mod_mat <- module_abundance(ps)
main_mod <- top_modules(mod_mat, 4)

mod_long <- as.data.frame(t(mod_mat)) |>
  rownames_to_column("sample_id") |>
  pivot_longer(-sample_id, names_to = "Module", values_to = "Abundance") |>
  left_join(meta |> select(sample_id, oceanic_sector, longitude), by = "sample_id") |>
  mutate(
    Module = ifelse(Module %in% main_mod, paste0("Module ", Module), "Other modules"),
    Module = factor(Module, levels = c(paste0("Module ", main_mod), "Other modules"))
  ) |>
  group_by(sample_id, oceanic_sector, longitude, Module) |>
  summarise(Abundance = sum(Abundance), .groups = "drop") |>
  arrange(longitude)

ggplot(mod_long, aes(longitude, Abundance, fill = Module)) +
  geom_area(color = "white", linewidth = 0.1) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  labs(x = "Longitude (°E)", y = "Module relative abundance") +
  theme_bw()

