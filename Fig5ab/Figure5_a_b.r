
rm(list=ls())

library(ggClusterNet)
library(phyloseq)
library(dplyr)
library(WGCNA)
library(sna)
library(igraph)
library(tidyverse)
library(network)
library(ggrepel)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(ggraph)



folder="16S_20251016"

ASV_all_16S_With_Tax <- read.table("ASV_Abundance_16S_split.txt", sep = "\t", header = TRUE, row.names = 1)
ASV_all_16S <- ASV_all_16S_With_Tax[, 1:(ncol(ASV_all_16S_With_Tax)-7)]
ASV_all_16S_Tax <- ASV_all_16S_With_Tax[ , (ncol(ASV_all_16S_With_Tax)-6):ncol(ASV_all_16S_With_Tax)]

infor <- data.frame(SampleID = colnames(ASV_all_16S), Group = NA, stringsAsFactors = FALSE)
infor$Group <- ifelse(substr(infor$SampleID,1,3) == "Bla", "Black",
                      ifelse(substr(infor$SampleID,1,3) == "Whi", "White",
                             ifelse(substr(infor$SampleID,1,3) == "Gre", "Green", "CK")))
row.names(infor) <- infor$SampleID


ASV_filtered_16S <- ASV_all_16S[rowSums(ASV_all_16S) > ncol(ASV_all_16S), ]
ASV_filtered_16S <- ASV_filtered_16S[, order(colnames(ASV_filtered_16S))]

common_rows_16S <- intersect(rownames(ASV_all_16S_Tax), rownames(ASV_filtered_16S))
ASV_filtered_16S_Tax <- ASV_all_16S_Tax[common_rows_16S, ]

ps_16S <- phyloseq(
  otu_table(as.matrix(ASV_filtered_16S), taxa_are_rows = TRUE),
  tax_table(as.matrix(ASV_filtered_16S_Tax)),
  sample_data(infor)
)



ps_all <- ps_16S
sample_data(ps_all)$Group <- factor("All")


set.seed(111)


tab.r_16S <- network.pip(
  ps = ps_all,
  big = TRUE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.75,
  p.threshold = 0.001,
  method = "spearman",
  label = TRUE,
  lab = "elements",
  group = "Group",
  fill = "Phylum",
  size = "igraph.degree",
  N = 400,
  ncpus = 8
)
saveRDS(tab.r_16S, file = paste0(folder,"/","16S_network_all.pip.rds"))

tab.r_16S<-readRDS(file = paste0(folder,"/","16S_network_all.pip.rds"))

dat_16S <- tab.r_16S[[2]]
cortab_16S <- dat_16S$net.cor.matrix$cortab
node_all <- dat_16S$net.cor.matrix$node
edge_all <- dat_16S$net.cor.matrix$edge



groups <- sort(unique(as.character(infor$Group))) # e.g. CK, Black, Green, White

present_ASVs <- lapply(groups, function(g){
  samp_ids <- rownames(infor)[infor$Group == g]
  asv_ids <- rownames(ASV_filtered_16S)[rowSums(ASV_filtered_16S[, samp_ids, drop = FALSE]) > 0]
  return(as.character(asv_ids))
})
names(present_ASVs) <- groups

nodes_list <- list()
edges_list <- list()
graphs_list <- list()
isolates_list <- list()



for (g in groups) {
  present <- present_ASVs[[g]]

  if ("OTU_1" %in% names(edge_all) & "OTU_2" %in% names(edge_all)) {
    edg <- edge_all %>% filter(OTU_1 %in% present & OTU_2 %in% present)
  } else if ("from" %in% names(edge_all) & "to" %in% names(edge_all)) {
    edg <- edge_all %>% filter(from %in% present & to %in% present)
  } else {
    edg <- cortab_16S %>% filter(OTU_1 %in% present & OTU_2 %in% present)
  }
  if ("elements" %in% names(node_all)) {
    nod <- node_all %>% filter(elements %in% present)
  } else if ("name" %in% names(node_all)) {
    nod <- node_all %>% filter(name %in% present)
  } else {
    stop("error")
  }

  


  if (nrow(edg) == 0) {
    warning(paste0("Group", g, " has no edge, skip"))
    nodes_list[[g]] <- nod
    edges_list[[g]] <- edg
    graphs_list[[g]] <- NULL
    next
  }
  
  if (all(c("OTU_1", "OTU_2") %in% names(edg))) {
    edg2 <- edg %>% rename(from = OTU_1, to = OTU_2)
  } else if (all(c("from", "to") %in% names(edg))) {
    edg2 <- edg
  } else if (all(c("source", "target") %in% names(edg))) {
    edg2 <- edg %>% rename(from = source, to = target)
  } else {
    stop(paste0("error"))
  }
  
  if ("elements" %in% names(nod)) {
    nod2 <- nod %>% rename(name = elements)
  } else if ("name" %in% names(nod)) {
    nod2 <- nod
  } else if ("id" %in% names(nod)) {
    nod2 <- nod %>% rename(name = id)
  } else {
    stop("error")
  }
  
  nod2$name <- as.character(nod2$name)
  edg2$from <- as.character(edg2$from)
  edg2$to   <- as.character(edg2$to)
  

  nod2<-nod2 %>% select(name,everything())
  
  edg2 <- edg2 %>% distinct(from, to, .keep_all = TRUE)
  edg2<-edg2 %>% select(from, to, everything())
  edg2 <- edg2 %>%
    filter((from %in% nod2$name) & (to %in% nod2$name))
  
  g_sub2 <- graph_from_data_frame(d = edg2, vertices = nod2, directed = FALSE)
  
  isolates <- V(g_sub2)[degree(g_sub2) == 0]
  if (length(isolates) > 0) {
    g_sub2 <- delete_vertices(g_sub2, isolates)
  }
  
  
  nodes_2 <- igraph::as_data_frame(g_sub2, what = "vertices")
  edges_2 <- igraph::as_data_frame(g_sub2, what = "edges") %>% 
    filter((from %in% nodes_2$name) & (to %in% nodes_2$name))
 
  nodes_3 <- nodes_2 %>% 
    select(name,Domain,Phylum,Class,
           Order,Family,Genus,Species)
  
  
  edges_3 <- edges_2 %>% 
    select(from, to, weight) %>% 
    mutate(cor = ifelse(weight>0, "Positive","Negative")) %>%
    mutate(weight = abs(weight))
  
  
  g_sub_3 <- graph_from_data_frame(d = edges_3, vertices = nodes_3, directed = FALSE)
  V(g_sub_3)$igraph.degree <- degree(g_sub_3)
  
  nodes_4 <- igraph::as_data_frame(g_sub_3, what = "vertices")
  edges_4 <- igraph::as_data_frame(g_sub_3, what = "edges")
  
  
  nodes_5 <- nodes_4 %>% 
    mutate(n_nodes = nrow(nodes_4)) %>% 
    mutate(n_edge = nrow(edges_4)) %>% 
    mutate(group_label = paste0(g, ": (nodes: ",n_nodes,"; links: ",n_edge,")" ))
  
  edges_5 <- edges_4 %>% 
    mutate(n_nodes = nrow(nodes_4)) %>% 
    mutate(n_edge = nrow(edges_4)) %>%
    mutate(group_label = paste0(g, ": (nodes: ",n_nodes,"; links: ",n_edge,")" ))
  
  g_sub_5 <- graph_from_data_frame(d = edges_5, vertices = nodes_5, directed = FALSE)
  
  
  nodes_list[[g]] <- nodes_5
  edges_list[[g]] <- edges_5
  graphs_list[[g]] <- g_sub_5
  
  isolates_list[[g]] <- isolates
  
}




View(graphs_list$Black)
View(nodes_list$Black)
View(edges_list$Black)
View(nodes_list$CK)
View(edges_list$CK)

vertex_attr_names(graphs_list$Black)
edge_attr_names(graphs_list$Black)
V(graphs_list$Black)$igraph.degree
E(graphs_list$Black)$cor


View(nodes_list$Black %>% select(name,Domain,Phylum,Class,
       Order,Family,Genus,Species)
)



isolates_df <- data.frame(
  group = names(isolates_list),
  n_isolates = sapply(isolates_list, length)
)





g<-graphs_list$Black

E(g)$cor


p_test <- ggraph(g, layout = "fr") +

  geom_edge_link(aes(colour = cor), alpha = 0.6, show.legend = TRUE) +
  
  geom_node_point(aes(fill = Phylum, size = igraph.degree), shape = 21, stroke = 0.3) +
  

scale_fill_manual(values = c("Proteobacteria"="#FB8072", "Actinobacteriota" ="#80B1D3",
                             "Cyanobacteria"="#FDB462", "Chloroflexi"="#8DD3C7",
                             "Acidobacteriota"="#FFFFB3","Bacteroidota"="#BEBADA",
                             "Planctomycetota"="#B3DE69", "Gemmatimonadota"="#FCCDE5",
                             "Myxococcota"="#BC80BD","Deinococcota"="#CCEBC5",
                             "Others"="grey80"),
                  name = "Phylum", na.value = "gray80") +


  scale_edge_colour_manual(values = c(Negative = "#6D98B5", Positive = "#D48852", zero = "gray70"),
                           name = "Correlation") +
  

  scale_size_continuous(range = c(2, 8), name = "Degree") +
  scale_edge_width(range = c(0.3, 1.5), guide = "none") +
  
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  ggtitle(V(g)$group_label[1])


p_test







phyla_colors<-c("Proteobacteria"="#FB8072", "Actinobacteriota" ="#80B1D3",
               "Cyanobacteria"="#FDB462", "Chloroflexi"="#8DD3C7",
               "Acidobacteriota"="#FFFFB3","Bacteroidota"="#BEBADA",
               "Planctomycetota"="#B3DE69", "Gemmatimonadota"="#FCCDE5",
               "Myxococcota"="#BC80BD","Deinococcota"="#CCEBC5",
               "Others"="grey80")

edge_colors <- c(Negative = "#6D98B5", Positive = "#D48852")


Samples_Palette <- c("CK"="#EEA236FF",
                     "Black"="#D43F3AFF",
                     "Green"="#5CB85CFF",
                     "White"="#46B8DAFF")

all_degrees <- unlist(lapply(graphs_list, function(g) V(g)$igraph.degree))
degree_range <- range(all_degrees)

size_breaks <- pretty(degree_range, n=5)
size_limits <- degree_range

desired_order <- c("CK", "Black", "Green", "White")


plot_list_unified <- list()

set.seed(123)
for (group_name in names(graphs_list)) {
  g <- graphs_list[[group_name]]
  
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(colour = cor), linewidth = 0.3, alpha = 0.6, show.legend = FALSE) +
    geom_node_point(aes(fill = Phylum, size = igraph.degree), shape = 21, stroke = 0.3, show.legend = FALSE) +
    scale_fill_manual(values = phyla_colors, name = "Phylum", na.value = "gray80") +
    scale_edge_colour_manual(values = edge_colors, name = "Correlation") +
    scale_size_continuous(breaks = size_breaks, limits = size_limits, range = c(1, 6), name = "Degree") +
    scale_edge_width(range = c(0.3, 1.5), guide = "none") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10)
    ) +
    ggtitle(V(g)$group_label[1])
  
  plot_list_unified[[group_name]] <- p
}

combined_plot_unified <- wrap_plots(
  plot_list_unified[desired_order],
  ncol = 4)

p_for_legend <- ggraph(graphs_list$Green, layout = "fr") +
  geom_edge_link(aes(colour = cor), linewidth = 0.3, alpha = 0.6, show.legend = TRUE) +
  geom_node_point(aes(fill = Phylum, size = igraph.degree), shape = 21, stroke = 0.3, show.legend = TRUE) +
  scale_fill_manual(values = phyla_colors, name = "Phylum", na.value = "gray80") +
  scale_edge_colour_manual(values = edge_colors, name = "Correlation") +
  scale_size_continuous(breaks = size_breaks, limits = size_limits, range = c(1, 6), name = "Degree") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.justification = "center",
        legend.box.just = "center"
) +
  guides(
    fill = guide_legend(nrow = 2, title.position = "left", title.hjust = 0.5),
    colour = guide_legend(nrow = 1, title.position = "left", title.hjust = 0.5),
    size = guide_legend(nrow = 1, title.position = "left", title.hjust = 0.5)
  )

legend_unified <- get_legend(p_for_legend)

final_plot_unified <- combined_plot_unified / legend_unified + 
  plot_layout(heights = c(10, 1))

print(final_plot_unified)

ggsave(
  filename = paste0(folder, "/16S_network_combined_unified_legend.pdf"),
  plot = final_plot_unified,
  width = 14,
  height = 4,
  dpi = 300
)








library(igraph)

quantify_network_looseness_full <- function(g, use_normalized_centralization = TRUE) {
  g <- simplify(as_undirected(g))
  
  n <- vcount(g)
  m <- ecount(g)
  
  is_conn <- is_connected(g)
  
  dens  <- edge_density(g)
  mdeg  <- if (n > 0) mean(degree(g)) else NA_real_
  clust <- transitivity(g, type = "average")
  cent  <- if (use_normalized_centralization) {
    centr_degree(g, normalized = TRUE)$centralization
  } else {
    centr_degree(g)$centralization
  }
  
  apl   <- if (is_conn && n > 1) mean_distance(g) else NA_real_
  diam  <- if (is_conn && n > 1) diameter(g)     else NA_real_
  
  comp <- components(g)
  main_ratio <- if (n > 0) max(comp$csize) / n else NA_real_
  
  neg_ratio <- sum(E(g)$cor == "Negative") / m
  
  data.frame(
    n_nodes              = n,
    n_edges              = m,
    edge_density         = dens,
    mean_degree          = mdeg,
    avg_path_length      = apl,
    diameter             = diam,
    clustering           = clust,
    degree_centralization= cent,
    main_component_ratio = main_ratio,
    is_connected         = is_conn,
    negative_link_ratio  = neg_ratio
    
  )
}

looseness_stats <- do.call(rbind, lapply(names(graphs_list), function(nm){
  out <- quantify_network_looseness_full(graphs_list[[nm]])
  cbind(Group = nm, out, row.names = NULL)
}))

print(looseness_stats)



looseness_stats_revised<-looseness_stats %>% 
  select(Group,
         clustering,
         degree_centralization,
         edge_density,
         mean_degree,
         negative_link_ratio
         ) %>% 
  mutate(negative_link_ratio = round(negative_link_ratio * 100, 2)) %>% 
  rename(`Negative correlation ratio (%)` = negative_link_ratio) %>% 
  mutate(clustering = round(clustering * 100, 2)) %>% 
  rename(`Clustering coefficient (%)` = clustering) %>% 
  mutate(edge_density = round(edge_density * 100, 2)) %>% 
  rename(`Edge density (%)` = edge_density) %>% 
  mutate(degree_centralization = round(degree_centralization * 100, 2)) %>% 
  rename(`Degree centralization (%)` = degree_centralization) %>% 
  mutate(mean_degree = round(mean_degree,2)) %>% 
  rename(`Mean degree`=mean_degree)
  
  




#### plot 
library(tidyverse)

looseness_long <- looseness_stats_revised %>%
  pivot_longer(
    cols = -Group,
    names_to = "Metric",
    values_to = "Value"
  )

looseness_long$Group <- factor(looseness_long$Group,
                               levels = c("CK", "Black", "Green", "White"))

p_all <- ggplot(looseness_long, aes(x = Group, y = Value, color = Group, fill = Group)) +
  geom_bar(stat = "identity", alpha = 0.2) +
  geom_text(aes(label = sprintf("%.2f", Value)), vjust = 1.5, size = 4) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 5) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  scale_color_manual(values = Samples_Palette) +
  scale_fill_manual(values = Samples_Palette) +
  labs(
    title = NULL,
    x = NULL,
    y = NULL
  )+
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    strip.placement = "outside"
  )

p_all

ggsave(paste0(folder, "/", "All_network_metrics_facet.pdf"),
       p_all, width = 14, height = 3)






library(ggpubr)

combined_ggpubr <- ggarrange(
  final_plot_unified, p_all,
  ncol = 1, nrow = 2,
  heights = c(4, 3)
)

combined_ggpubr

ggsave(file.path(folder, "combined_ggpubr_16S.pdf"),
       combined_ggpubr, width = 14, height = 7)





