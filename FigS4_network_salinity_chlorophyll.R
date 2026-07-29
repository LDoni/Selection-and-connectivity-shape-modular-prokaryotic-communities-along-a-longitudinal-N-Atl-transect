## Figure S4 - network salinity and chlorophyll

source("functions/01_atlantic_transect_and_community.R")
source("functions/03_network_and_modules.R")

g <- read_graph("output/raw_tables/SparCC_Network_pos.graphml", format = "graphml")
lay <- network_layout(g)
v_size <- network_node_size(g)


### Fig S4A Salinity

plot(g, layout = lay, vertex.size = v_size, vertex.label = NA,
     vertex.color = network_continuous_colors(V(g)$salinity_psu, "RdYlBu", reverse = TRUE),
     edge.color = adjustcolor("grey70", alpha.f = 0.25), main = "Salinity")
network_continuous_scale(V(g)$salinity_psu, "RdYlBu", reverse = TRUE, label = "PSU")


### Fig S4B chl

plot(g, layout = lay, vertex.size = v_size, vertex.label = NA,
     vertex.color = network_continuous_colors(V(g)$chlorophyll_mg_m3, "Greens"),
     edge.color = adjustcolor("grey70", alpha.f = 0.25), main = "Chlorophyll")
network_continuous_scale(V(g)$chlorophyll_mg_m3, "Greens", label = expression("mg m"^-3))


