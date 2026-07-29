## SparCC network and modules

library(dplyr)
library(tibble)
library(igraph)
source("functions/01_atlantic_transect_and_community.R")

cor_file <- "data/network/Cor_SparCC_Prok_all.csv"
pval_file <- "data/network/Pval_SparCC_Prok_all.csv"


  ##  SparCC matrices are not included in the repository. 
  #  Put Cor_SparCC_Prok_all.csv and Pval_SparCC_Prok_all.csv in data/network/ 
   # to rebuild output/raw_tables/SparCC_Network_pos.graphml and output/raw_tables/cluster.rds


cor_mat <- as.matrix(read.csv(cor_file, row.names = 1, check.names = FALSE))
p_mat <- as.matrix(read.csv(pval_file, row.names = 1, check.names = FALSE))

idx <- which(upper.tri(cor_mat), arr.ind = TRUE)
edges <- data.frame(
  from = rownames(cor_mat)[idx[, 1]],
  to = colnames(cor_mat)[idx[, 2]],
  r = cor_mat[idx],
  p = p_mat[idx]
) |>
  mutate(fdr = p.adjust(p, method = "BH")) |>
  filter(r >= 0.51, fdr < 0.05)

g <- graph_from_data_frame(edges, directed = FALSE)
cl <- cluster_edge_betweenness(g)

V(g)$Cluster <- as.character(membership(cl)[V(g)$name])
V(g)$Degree <- degree(g)
V(g)$n <- as.numeric(sizes(cl)[membership(cl)[V(g)$name]])

ps_raw <- read_atlantic_raw_ps()
ps_rel <- tss_normalize_ps(ps_raw)
otu_rel <- as(otu_table(ps_rel), "matrix")
otu_raw <- as(otu_table(ps_raw), "matrix")
if (!taxa_are_rows(ps_rel)) otu_rel <- t(otu_rel)
if (!taxa_are_rows(ps_raw)) otu_raw <- t(otu_raw)

otu_rel <- otu_rel[V(g)$name, , drop = FALSE]
otu_raw <- otu_raw[V(g)$name, , drop = FALSE]
max_sample_index <- max.col(otu_rel, ties.method = "first")
max_sample <- colnames(otu_rel)[max_sample_index]
meta <- sample_metadata(ps_rel)
node_meta <- meta[match(max_sample, meta$sample_id), ]

V(g)$max_relative_abundance <- otu_rel[cbind(seq_len(nrow(otu_rel)), max_sample_index)]
V(g)$total_counts <- rowSums(otu_raw)
V(g)$sample_id <- max_sample
for (v in c("day_of_navigation", "sampling_date", "oceanic_sector", "latitude", "longitude",
            "sst_degC", "salinity_psu", "chlorophyll_mg_m3", "oxygen_mmol_m3",
            "nitrate_mmol_m3", "phosphate_mmol_m3", "silicate_mmol_m3")) {
  value <- node_meta[[v]]
  if (inherits(value, "Date") || is.factor(value)) value <- as.character(value)
  g <- set_vertex_attr(g, v, value = value)
}

cluster_table <- data.frame(
  Cluster = V(g)$Cluster,
  OTU_ID = V(g)$name,
  n = V(g)$n,
  Degree = V(g)$Degree
) |>
  arrange(as.numeric(Cluster), OTU_ID)

write_graph(g, "output/raw_tables/SparCC_Network_pos.graphml", format = "graphml")
saveRDS(cluster_table, "output/raw_tables/cluster.rds")
write.csv(edges, "output/raw_tables/network_edges_filtered.csv", row.names = FALSE)

data.frame(
  nodes = vcount(g),
  edges = ecount(g),
  modularity = modularity(cl),
  modules = length(unique(cluster_table$Cluster))
)
