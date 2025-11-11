######Ferroptosis drives an immunosuppressive microenvironment through reprogramming
###################################################################


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "msigdbr", "enrichplot", "ggplot2", "dplyr","readxl",
                       "org.Hs.eg.db","STRINGdb","igraph","ggraph","fmsb","tibble","tidyr","readr"))

library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(readxl)
library(org.Hs.eg.db)
library(STRINGdb)
library(igraph)
library(ggraph)
library(fmsb)
library(tibble)
library(tidyr)
library(readr)



output_dir <- "my_results/immune"



#### Ferroptosis Mechanism of High Activity Cells Macro and Its Association with Immune Function
###GSEA identifies key pathways
diff_genes <- read_excel("my_data/diff_genes_PAAD_immune_detailed_Macro.xlsx")

gene_rank <- diff_genes$logFC
names(gene_rank) <- diff_genes$Gene
gene_rank <- sort(gene_rank, decreasing = TRUE)

all_sets <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  filter(grepl("FERR|IRON|GLUTATHIONE|LIPID|FATTY_ACID|PEROXID|OXIDATIVE|GSH|INFLAMM|IMMUN|CYTOKINE|CHEMOKINE|INTERLEUKIN|INTERFERON", 
               gs_name, ignore.case = TRUE))
set.seed(123) 

#GSEA analysis of ferroptosis-related pathways
ferro_sets <- all_sets %>% 
  filter(grepl("FERR|IRON|GLUTATHIONE|LIPID|FATTY_ACID|PEROXID|OXIDATIVE|GSH", 
               gs_name, ignore.case = TRUE))
gsea_ferro <- GSEA(gene_rank, 
                   TERM2GENE = ferro_sets[, c("gs_name", "gene_symbol")],
                   pvalueCutoff = 0.05)
write.csv(as.data.frame(gsea_ferro), file= file.path(output_dir, "Macro_GSEA_Ferroptosisresults.csv"), row.names = FALSE)
pdf(file.path(output_dir, "Macro_GSEA_FerroptosisPathways.pdf"), width = 12, height = 8)
gseaplot2(gsea_ferro, 
          geneSetID = 1:min(3, nrow(gsea_ferro)),
          pvalue_table = TRUE,
          title = "Ferroptosis-related Pathways in Macro Cells",
          color = c("#E64B35", "#4DBBD5", "#00A087"),
          base_size = 14,
          ES_geom = "line")
dev.off()

#GSEA analysis of immune-related pathways
immune_sets <- all_sets %>% 
  filter(grepl("INFLAMM|IMMUN|CYTOKINE|CHEMOKINE|INTERLEUKIN|INTERFERON", 
               gs_name, ignore.case = TRUE))
gsea_immune <- GSEA(gene_rank, 
                    TERM2GENE = immune_sets[, c("gs_name", "gene_symbol")],
                    pvalueCutoff = 0.05)
write.csv(as.data.frame(gsea_immune), file= file.path(output_dir, "Macro_GSEA_Immuneresults.csv"), row.names = FALSE)
pdf(file.path(output_dir, "Macro_GSEA_ImmunePathways.pdf"), width = 12, height = 8)
gseaplot2(gsea_immune, 
          geneSetID = 1:min(3, nrow(gsea_immune)),
          pvalue_table = TRUE,
          title = "Immune-related Pathways in Macro Cells",
          color = c("#3C5488", "#F39B7F", "#8491B4"),
          base_size = 14,
          ES_geom = "line")
dev.off()



###Network analysis to identify key driver genes
gsea_results <- read.csv("my_data/Macro_GSEA_results.csv")
sig_pathways <- subset(gsea_results, p.adjust < 0.05 & abs(NES) > 1.5)

core_genes <- unique(unlist(lapply(strsplit(sig_pathways$core_enrichment, "/"), function(x) {
  x[1:round(length(x)*0.2)]        
})))

ferr_db_genes <- read.csv("my_data/ferroptosis_all.csv")
known_ferro_genes <- unique(ferr_db_genes$symbol)

#Screen candidate driver genes
#Criteria: P-value < 0.05, |logFC| > 1, and present in the core gene set or known ferroptosis gene set
candidate_driver_genes <- diff_genes %>% 
  filter(PValue < 0.05 & abs(logFC) > 1) %>% 
  filter(Gene %in% c(core_genes, known_ferro_genes)) %>% 
  arrange(desc(abs(logFC)))

#Construct PPI network
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400,input_directory="", protocol="http")

driver_genes_mapped <- string_db$map(
  data.frame(Gene = candidate_driver_genes$Gene),  
  "Gene", 
  removeUnmappedRows = TRUE
)

interactions <- string_db$get_interactions(driver_genes_mapped$STRING_id) %>%
  distinct(from, to, .keep_all = TRUE) %>%
  mutate(
    from = gsub("^9606\\.", "", from),
    to = gsub("^9606\\.", "", to)
  )

driver_genes_mapped <- driver_genes_mapped %>%
  mutate(STRING_id = gsub("^9606\\.", "", STRING_id))

valid_interactions <- interactions %>%
  filter(from %in% driver_genes_mapped$STRING_id & 
           to %in% driver_genes_mapped$STRING_id)

net <- graph_from_data_frame(
  d = valid_interactions[, c("from", "to", "combined_score")],
  directed = FALSE
)

V(net)$Gene <- driver_genes_mapped$Gene[match(V(net)$name, driver_genes_mapped$STRING_id)]

#Compute Node Importance Metric
V(net)$degree <- degree(net)
V(net)$betweenness <- betweenness(net)
V(net)$hub_score <- hub_score(net)$vector
V(net)$closeness <- closeness(net)
V(net)$eigen_centrality <- eigen_centrality(net)$vector
V(net)$page_rank <- page_rank(net)$vector

#Identify key driver genes using PCA
vertex_data <- data.frame(
  STRING_id = V(net)$name,
  Gene = V(net)$Gene,
  degree = V(net)$degree,
  betweenness = V(net)$betweenness,
  closeness = V(net)$closeness,
  eigen_centrality = V(net)$eigen_centrality,
  page_rank = V(net)$page_rank,
  hub_score = V(net)$hub_score
) 

scale_metrics <- function(x) {
  (x - mean(x)) / sd(x)
}

vertex_data_norm <- vertex_data %>%
  mutate(across(c(degree:hub_score), scale_metrics))

pca_res <- prcomp(vertex_data_norm[, c("degree", "betweenness", "closeness", "eigen_centrality", "page_rank","hub_score")], 
                  scale. = FALSE)
vertex_data_norm$pca_score <- pca_res$x[,1]  

#Save the list of key driver genes
hub_genes_pca <- vertex_data_norm %>%
  arrange(desc(pca_score))

write.csv(hub_genes_pca, file= file.path(output_dir, "Macro_Ferroptosis_DriverGenes_Detailed.csv"), row.names = FALSE)

hub_genes <- vertex_data_norm %>%
  arrange(desc(pca_score)) %>%
  head(10) %>%
  pull(Gene)

driver_genes <- candidate_driver_genes %>% 
  mutate(
    is_hub = Gene %in% hub_genes,
    in_core_pathway = Gene %in% core_genes,
    in_known_ferro = Gene %in% known_ferro_genes
  )

write.csv(driver_genes, file= file.path(output_dir, "Macro_Ferroptosis_DriverGenes_Type.csv"), row.names = FALSE)

##Visual Network
pdf(file.path(output_dir, "Macro_DriverGenes_PPI_Network.pdf"), width = 14, height = 12)
set.seed(123)
ggraph(net, layout = "kk") + 
  geom_edge_diagonal(aes(width = combined_score, 
                         alpha = combined_score),
                     colour = "#666666",
                     strength = 0.8,
                     show.legend = FALSE) + 
  geom_node_point(aes(size = degree, 
                      color = ifelse(Gene %in% hub_genes, "Hub Gene", "Other Gene"),
                      fill = vertex_data_norm$pca_score), 
                  shape = 21, stroke = 1.5, alpha = 1) +
  geom_node_text(aes(label = Gene, 
                     filter = Gene %in% hub_genes),  
                 size = 10, color = "#e12216", fontface = "bold",
                 repel = TRUE, box.padding = 0.3, max.overlaps = 100) +
  geom_node_text(aes(label = Gene, 
                     filter = !Gene %in% hub_genes),  
                 size = 10, color = "black", fontface = "bold",
                 repel = TRUE, box.padding = 0.3, max.overlaps = 100) +
  scale_edge_width(range = c(3, 8)) +
  scale_size_continuous(range = c(4, 15),
                        breaks = seq(min(V(net)$degree), max(V(net)$degree), length.out = 4),
                        name = "Degree") +
  scale_fill_gradientn(colors = c("#f7f100", "#e12216"),
                       name = "PCA Score") +
  scale_color_manual(values = c("Hub Gene" = "#e12216", "Other Gene" = "black"),
                     name = "Gene Type") +
  guides(color = guide_legend(override.aes = list(size = 5)),
         size = guide_legend(),
         fill = guide_colorbar()) +
  labs(title = "Macrophage Ferroptosis Driver Genes PPI Network",
       subtitle = paste("Top 10 Driver Genes | Ferroptosis Score =", round(ferro_score, 3))) +
  theme_void() +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        plot.caption = element_text(size = 10, hjust = 0.5, margin = margin(t = 10)),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.spacing.x = unit(0.5, 'cm'))
dev.off()

#Radar chart visualizing key driver genes across multiple dimensions
radar_data <- vertex_data_norm %>%
  dplyr::filter(Gene %in% hub_genes) %>%
  dplyr::select(Gene, degree, betweenness, closeness, eigen_centrality, page_rank, hub_score) %>%
  tibble::column_to_rownames("Gene") %>%
  as.data.frame()

radar_data <- rbind(rep(1,6), rep(0,6), radar_data)

pdf(file.path(output_dir, "Macro_DriverGenes_RadarPlot.pdf"), width = 10, height = 10)
par(mar = c(1,1,3,1))
layout(matrix(1:4, nrow = 2, byrow = TRUE))
colors <- rainbow(length(hub_genes))
for(i in 1:length(hub_genes)){
  gene <- hub_genes[i]
  data <- radar_data[c(1,2, which(rownames(radar_data) == gene)), ]
  
  radarchart(data,
             axistype = 1,
             pcol = colors[i], pfcol = paste0(colors[i], "60"), plwd = 2,
             cglcol = "grey", cglty = 1, cglwd = 0.8,
             axislabcol = "grey20",
             vlcex = 0.8, title = gene)
}
dev.off()



###Predict TF
top_drivers <- hub_genes
trrust_raw <- read_tsv("my_data/trrust_rawdata.human.tsv", col_names = c("TF", "Target", "Direction", "PMID"))

#Screen TFs related to key driver genes
tf_network <- trrust_raw %>% 
  filter(Target %in% top_drivers) %>% 
  distinct(TF, Target, .keep_all = TRUE)

#Construct TF-target gene networks
tf_graph <- graph_from_data_frame(
  tf_network[, c("TF", "Target")],
  directed = TRUE
)

# Calculate network centrality metrics
V(tf_graph)$type <- ifelse(V(tf_graph)$name %in% tf_network$TF, "TF", "Target")
V(tf_graph)$degree <- degree(tf_graph)
V(tf_graph)$betweenness <- betweenness(tf_graph)
V(tf_graph)$closeness <- closeness(tf_graph)
V(tf_graph)$eigen_centrality <- eigen_centrality(tf_graph)$vector
V(tf_graph)$page_rank <- page_rank(tf_graph)$vector
V(tf_graph)$hub_score <- hub_score(tf_graph)$vector

#Identify key TFs using PCA
vertex_data <- data.frame(
  Gene = V(tf_graph)$name,
  type = V(tf_graph)$type,
  degree = V(tf_graph)$degree,
  closeness = V(tf_graph)$closeness,
  eigen_centrality = V(tf_graph)$eigen_centrality,
  hub_score = V(tf_graph)$hub_score
) 

tf_data <- vertex_data %>% 
  filter(type == "TF")

scale_metrics <- function(x) {
  (x - mean(x)) / sd(x)
}
tf_data_norm <- tf_data %>%
  mutate(across(c(degree:hub_score), scale_metrics))

pca_res <- prcomp(tf_data_norm[, c("degree", "closeness", 
                                   "eigen_centrality",  "hub_score")], 
                  scale. = FALSE)
tf_data_norm$pca_score <- pca_res$x[,1]

#Save key TF information
tf_data_norm <- tf_data_norm %>%
  arrange(desc(pca_score))

write.csv(tf_data_norm, file= file.path(output_dir, "Macro_Ferroptosis_KeyTFs.csv"), row.names = FALSE)

top_tfs <- tf_data_norm %>%
  head(10) %>%
  pull(Gene)

#Visualize TF-target gene network
sub_nodes <- unique(c(top_tfs, top_drivers))
sub_edges <- tf_network %>% 
  filter(TF %in% sub_nodes & Target %in% sub_nodes)
sub_graph <- graph_from_data_frame(sub_edges, directed = TRUE)

V(sub_graph)$pca_score <- tf_data_norm$pca_score[match(V(sub_graph)$name, tf_data_norm$Gene)]
V(sub_graph)$node_type <- ifelse(V(sub_graph)$name %in% top_tfs, "Key TF", 
                                 ifelse(V(sub_graph)$name %in% top_drivers, "Driver Gene", "Other"))
V(sub_graph)$size <- ifelse(V(sub_graph)$name %in% top_tfs, 
                            scales::rescale(V(sub_graph)$pca_score, to = c(8, 15)), 
                            8)
V(sub_graph)$color <- case_when(
  V(sub_graph)$node_type == "Key TF" ~ "#154996",
  V(sub_graph)$node_type == "Driver Gene" ~ "#e12216",
  TRUE ~ "#91D1C2"
)

pdf(file.path(output_dir, "Macro_Ferroptosis_TF_Network.pdf"), width = 14, height = 12)
set.seed(123)
ggraph(sub_graph, layout = "kk") + 
  geom_edge_link(aes(color = "#666666"), 
                 arrow = arrow(length = unit(5, 'mm')), 
                 end_cap = circle(5, 'mm'),
                 alpha = 1) +
  geom_node_point(aes(size = size, color = color, alpha = ifelse(node_type == "Key TF", pca_score, 1)), 
                  show.legend = FALSE) +
  geom_node_text(aes(label = name, 
                     filter = node_type %in% c("Key TF", "Driver Gene")),
                 size = 10, 
                 repel = FALSE,  
                 nudge_y = 0.1,  
                 nudge_x = 0.1,  
                 check_overlap = TRUE,  
                 box.padding = 0) + 
  scale_size_identity() +
  scale_color_identity() +
  scale_alpha_continuous(range = c(0.5, 1)) +
  labs(title = "Macrophage Ferroptosis: Key Transcription Factors and Driver Genes") +
  theme_void() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5))
dev.off()





#### Ferroptosis Mechanism of High Activity Cells CD8Tcm and Its Association with Immune Function
##GSEA identifies key pathways
diff_genes_CD8Tcm <- read_excel("my_data/diff_genes_PAAD_immune_detailed_CD8Tcm.xlsx")

gene_rank_CD8Tcm <- diff_genes_CD8Tcm$logFC
names(gene_rank_CD8Tcm) <- diff_genes_CD8Tcm$Gene
gene_rank_CD8Tcm <- sort(gene_rank_CD8Tcm, decreasing = TRUE)

all_sets <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  filter(grepl("FERR|IRON|GLUTATHIONE|LIPID|FATTY_ACID|PEROXID|OXIDATIVE|GSH|INFLAMM|IMMUN|CYTOKINE|CHEMOKINE|INTERLEUKIN|INTERFERON", 
               gs_name, ignore.case = TRUE))
set.seed(123) 

#GSEA analysis of ferroptosis-related pathways
ferro_sets <- all_sets %>% 
  filter(grepl("FERR|IRON|GLUTATHIONE|LIPID|FATTY_ACID|PEROXID|OXIDATIVE|GSH", 
               gs_name, ignore.case = TRUE))
gsea_ferro_CD8Tcm <- GSEA(gene_rank_CD8Tcm, 
                          TERM2GENE = ferro_sets[, c("gs_name", "gene_symbol")],
                          pvalueCutoff = 0.05)
write.csv(as.data.frame(gsea_ferro_CD8Tcm), file= file.path(output_dir, "CD8Tcm_GSEA_Ferroptosisresults.csv"), row.names = FALSE)
pdf(file.path(output_dir, "CD8Tcm_GSEA_FerroptosisPathways.pdf"), width = 12, height = 8)
gseaplot2(gsea_ferro_CD8Tcm, 
          geneSetID = 1:min(3, nrow(gsea_ferro_CD8Tcm)),
          pvalue_table = TRUE,
          title = "Ferroptosis-related Pathways in CD8Tcm Cells",
          color = c("#E64B35", "#4DBBD5", "#00A087"),
          base_size = 14,
          ES_geom = "line")
dev.off()

#GSEA analysis of immune-related pathways
immune_sets <- all_sets %>% 
  filter(grepl("INFLAMM|IMMUN|CYTOKINE|CHEMOKINE|INTERLEUKIN|INTERFERON", 
               gs_name, ignore.case = TRUE))
gsea_immune_CD8Tcm <- GSEA(gene_rank_CD8Tcm, 
                           TERM2GENE = immune_sets[, c("gs_name", "gene_symbol")],
                           pvalueCutoff = 0.05)
write.csv(as.data.frame(gsea_immune_CD8Tcm), file= file.path(output_dir, "CD8Tcm_GSEA_Immuneresults.csv"), row.names = FALSE)
pdf(file.path(output_dir, "CD8Tcm_GSEA_ImmunePathways.pdf"), width = 12, height = 8)
gseaplot2(gsea_immune_CD8Tcm, 
          geneSetID = 1:min(3, nrow(gsea_immune_CD8Tcm)),
          pvalue_table = TRUE,
          title = "Immune-related Pathways in CD8Tcm Cells",
          color = c("#3C5488", "#F39B7F", "#8491B4"),
          base_size = 14,
          ES_geom = "line")
dev.off()



###Network analysis to identify key driver genes
gsea_results_CD8Tcm <- read.csv("my_data/CD8Tcm_GSEA_results.csv")
sig_pathways_CD8Tcm <- subset(gsea_results_CD8Tcm, p.adjust < 0.05)

core_genes_CD8Tcm <- unique(unlist(lapply(strsplit(sig_pathways_CD8Tcm$core_enrichment, "/"), function(x) {
  x[1:round(length(x)*0.2)]        
})))

#Screen candidate driver genes
#Criteria: P-value < 0.05, |logFC| > 1, and present in the core gene set or known ferroptosis gene set
candidate_driver_genes_CD8Tcm <- diff_genes_CD8Tcm %>% 
  filter(PValue < 0.05 & abs(logFC) > 1) %>% 
  filter(Gene %in% c(core_genes_CD8Tcm, known_ferro_genes)) %>% 
  arrange(desc(abs(logFC)))

#Construct PPI network
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400,input_directory="", protocol="http")

driver_genes_mapped_CD8Tcm <- string_db$map(
  data.frame(Gene = candidate_driver_genes_CD8Tcm$Gene),  
  "Gene", 
  removeUnmappedRows = TRUE
)

interactions_CD8Tcm <- string_db$get_interactions(driver_genes_mapped_CD8Tcm$STRING_id) %>%
  distinct(from, to, .keep_all = TRUE) %>%
  mutate(
    from = gsub("^9606\\.", "", from),
    to = gsub("^9606\\.", "", to)
  )

driver_genes_mapped_CD8Tcm <- driver_genes_mapped_CD8Tcm %>%
  mutate(STRING_id = gsub("^9606\\.", "", STRING_id))

valid_interactions_CD8Tcm <- interactions_CD8Tcm %>%
  filter(from %in% driver_genes_mapped_CD8Tcm$STRING_id & 
           to %in% driver_genes_mapped_CD8Tcm$STRING_id)

#(No valid edges were found, so there is no need to construct a PPI network. The driver gene is the hub gene.)

driver_genes_CD8Tcm <- candidate_driver_genes_CD8Tcm %>% 
  mutate(
    in_core_pathway = Gene %in% core_genes_CD8Tcm,
    in_known_ferro = Gene %in% known_ferro_genes
  )
write.csv(driver_genes_CD8Tcm, file= file.path(output_dir, "CD8Tcm_Ferroptosis_DriverGenes_Type.csv"), row.names = FALSE)



###Predict TF
top_drivers_CD8Tcm <- driver_genes_CD8Tcm$Gene

#Screen TFs related to key driver genes
tf_network_CD8Tcm <- trrust_raw %>% 
  filter(Target %in% top_drivers_CD8Tcm) %>% 
  distinct(TF, Target, .keep_all = TRUE)

#Construct TF-target gene networks
tf_graph_CD8Tcm <- graph_from_data_frame(
  tf_network_CD8Tcm[, c("TF", "Target")],
  directed = TRUE
)

#Calculate network centrality metrics
V(tf_graph_CD8Tcm)$type <- ifelse(V(tf_graph_CD8Tcm)$name %in% tf_network_CD8Tcm$TF, "TF", "Target")
V(tf_graph_CD8Tcm)$degree <- degree(tf_graph_CD8Tcm)
V(tf_graph_CD8Tcm)$betweenness <- betweenness(tf_graph_CD8Tcm)
V(tf_graph_CD8Tcm)$closeness <- closeness(tf_graph_CD8Tcm)
V(tf_graph_CD8Tcm)$eigen_centrality <- eigen_centrality(tf_graph_CD8Tcm)$vector
V(tf_graph_CD8Tcm)$page_rank <- page_rank(tf_graph_CD8Tcm)$vector
V(tf_graph_CD8Tcm)$hub_score <- hub_score(tf_graph_CD8Tcm)$vector

#Identify key TFs using PCA
vertex_data_CD8Tcm <- data.frame(
  Gene = V(tf_graph_CD8Tcm)$name,
  type = V(tf_graph_CD8Tcm)$type,
  degree = V(tf_graph_CD8Tcm)$degree,
  betweenness = V(tf_graph_CD8Tcm)$betweenness,
  closeness = V(tf_graph_CD8Tcm)$closeness,
  eigen_centrality = V(tf_graph_CD8Tcm)$eigen_centrality,
  page_rank = V(tf_graph_CD8Tcm)$page_rank,
  hub_score = V(tf_graph_CD8Tcm)$hub_score
) 

#(There is only one TF and one target gene, so there is no need to use PCA to select high-scoring TFs)
#Save key TF information
tf_data_CD8Tcm <- vertex_data_CD8Tcm %>% 
  filter(type == "TF")
write.csv(tf_data_CD8Tcm, file= file.path(output_dir, "CD8Tcm_Ferroptosis_KeyTFs.csv"), row.names = FALSE)

top_tfs_CD8Tcm <- tf_data_CD8Tcm %>%
  head(10) %>%
  pull(Gene)

#Visualize TF-target gene network
sub_nodes_CD8Tcm <- unique(c(top_tfs_CD8Tcm, top_drivers_CD8Tcm))
sub_edges_CD8Tcm <- tf_network_CD8Tcm %>% 
  filter(TF %in% sub_nodes_CD8Tcm & Target %in% sub_nodes_CD8Tcm)
sub_graph_CD8Tcm <- graph_from_data_frame(sub_edges_CD8Tcm, directed = TRUE)

V(sub_graph_CD8Tcm)$node_type <- ifelse(V(sub_graph_CD8Tcm)$name %in% top_tfs_CD8Tcm, "Key TF", 
                                        ifelse(V(sub_graph_CD8Tcm)$name %in% top_drivers_CD8Tcm, "Driver Gene", "Other"))
V(sub_graph_CD8Tcm)$size <- 8
V(sub_graph_CD8Tcm)$color <- case_when(
  V(sub_graph_CD8Tcm)$node_type == "Key TF" ~ "#154996",
  V(sub_graph_CD8Tcm)$node_type == "Driver Gene" ~ "#e12216",
  TRUE ~ "#91D1C2"
)

pdf(file.path(output_dir, "CD8Tcm_Ferroptosis_TF_Network.pdf"), width = 14, height = 12)
set.seed(123)
ggraph(sub_graph_CD8Tcm, layout = "kk") + 
  geom_edge_link(aes(color = "#666666"), 
                 arrow = arrow(length = unit(5, 'mm')), 
                 end_cap = circle(5, 'mm'),
                 alpha = 1) +
  geom_node_point(aes(size = size, color = color, alpha = 1), 
                  show.legend = FALSE) +
  geom_node_text(aes(label = name, 
                     filter = node_type %in% c("Key TF", "Driver Gene")),
                 size = 10, 
                 repel = FALSE,  
                 nudge_y = 0.1,  
                 nudge_x = 0.1,  
                 check_overlap = TRUE,  
                 box.padding = 0) + 
  scale_size_identity() +
  scale_color_identity() +
  scale_alpha_continuous(range = c(0.5, 1)) +
  labs(title = "CD8Tcm Ferroptosis: Key Transcription Factors and Driver Genes") +
  theme_void() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5))

dev.off()





