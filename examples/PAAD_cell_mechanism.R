######Ferroptosis mechanisms in different cell types
###################################################################



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "msigdbr", "enrichplot", "ggplot2", "dplyr","readxl",
                       "DESeq2","org.Hs.eg.db","ggrepel","pheatmap","VennDiagram","grid",
                       "openxlsx","STRINGdb","igraph","ggraph","fmsb","tibble","tidyr","readr"))

library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(readxl)
library(DESeq2)
library(org.Hs.eg.db)
library(ggrepel)
library(pheatmap)
library(VennDiagram)
library(grid)
library(openxlsx)
library(STRINGdb)
library(igraph)
library(ggraph)
library(fmsb)
library(tibble)
library(tidyr)
library(readr)



output_dir <- "my_results/cellmechanism"



#### Analysis of the Functional Regulatory Mechanisms of Cells: cDC1 vs. cDC2
## Identify DEGs
# Load and Filter Differential Gene Data
cDC1_diff <- read_excel("my_data/diff_genes_PAAD_immune_detailed_cDC1.xlsx")
cDC2_diff <- read_excel("my_data/diff_genes_PAAD_immune_detailed_cDC2.xlsx")

cDC1_diff$cell_type <- "cDC1"
cDC2_diff$cell_type <- "cDC2"

combined_diff <- bind_rows(cDC1_diff, cDC2_diff)
sig_genes <- combined_diff %>% 
  filter(PValue < 0.05 & abs(logFC) > 1)

#Obtain the specific differential genes of the two cell types
cDC1_specific <- setdiff(
  sig_genes$Gene[sig_genes$cell_type == "cDC1"],
  sig_genes$Gene[sig_genes$cell_type == "cDC2"]
)

cDC2_specific <- setdiff(
  sig_genes$Gene[sig_genes$cell_type == "cDC2"],
  sig_genes$Gene[sig_genes$cell_type == "cDC1"]
)

#Veen plot
pdf(file.path(output_dir, "cDC2_diffgene_venn.pdf"), width = 14, height = 12)
venn.plot <- venn.diagram(
  x = list(
    cDC1 = sig_genes$Gene[sig_genes$cell_type == "cDC1"],
    cDC2 = sig_genes$Gene[sig_genes$cell_type == "cDC2"]
  ),
  filename = NULL,
  fill = c("#1b9e77", "#d95f02"),
  main = "Differential Genes Venn Diagram",
  main.cex = 2,
  cat.cex = 1.5,
  cex = 1.5
)
grid.draw(venn.plot)
dev.off()


#Output the specific differential genes of the two cell types and the information of common differential genes
common_genes <- intersect(
  sig_genes$Gene[sig_genes$cell_type == "cDC1"],
  sig_genes$Gene[sig_genes$cell_type == "cDC2"]
)

create_sorted_df <- function(genes, cell_type) {
  df <- sig_genes %>% 
    filter(Gene %in% genes, cell_type == !!cell_type) %>%
    arrange(desc(abs(logFC))) %>%  
    select(Gene, logFC, PValue) %>%
    mutate(importance_rank = row_number())  
  
  return(df)
}

cDC1_specific_df <- create_sorted_df(cDC1_specific, "cDC1")
cDC2_specific_df <- create_sorted_df(cDC2_specific, "cDC2")

common_df_cDC1 <- create_sorted_df(common_genes, "cDC1") %>%
  rename(logFC_cDC1 = logFC, PValue_cDC1 = PValue)
common_df_cDC2 <- create_sorted_df(common_genes, "cDC2") %>%
  rename(logFC_cDC2 = logFC, PValue_cDC2 = PValue)

common_df <- common_df_cDC1 %>%
  inner_join(common_df_cDC2, by = "Gene") %>%
  mutate(logFC_diff = abs(logFC_cDC1) - abs(logFC_cDC2)) %>%
  arrange(desc(abs(logFC_cDC1) + abs(logFC_cDC2)))  # Sort by the sum of the absolute values of logFC

wb <- createWorkbook()
addWorksheet(wb, "cDC1_specific")
addWorksheet(wb, "cDC2_specific")
addWorksheet(wb, "Common_genes")
writeData(wb, "cDC1_specific", cDC1_specific_df)
writeData(wb, "cDC2_specific", cDC2_specific_df)
writeData(wb, "Common_genes", common_df)

saveWorkbook(wb, file = file.path(output_dir, "cDC1&cDC2_diffgene_sc.xlsx"), overwrite = TRUE)


#Obtain ferroptosis-related genes among the specific differential genes of the two cell types
ferr_db_genes <- read.csv("my_data/ferroptosis_all.csv")
ferro_genes <- unique(ferr_db_genes$symbol)

cDC1_ferro_genes <- intersect(cDC1_specific, ferro_genes)
cDC1_ferro_ratio <- length(cDC1_ferro_genes) / length(cDC1_specific)

cDC2_ferro_genes <- intersect(cDC2_specific, ferro_genes)
cDC2_ferro_ratio <- length(cDC2_ferro_genes) / length(cDC2_specific)

ferro_results <- data.frame(
  CellType = c("cDC1", "cDC2"),
  Total_Specific_Genes = c(length(cDC1_specific), length(cDC2_specific)),
  Ferroptosis_Genes = c(length(cDC1_ferro_genes), length(cDC2_ferro_genes)),
  Ferroptosis_Ratio = c(cDC1_ferro_ratio, cDC2_ferro_ratio),
  Ferroptosis_Gene_List = c(paste(cDC1_ferro_genes, collapse = "; "), 
                            paste(cDC2_ferro_genes, collapse = "; "))
)
write.csv(ferro_results, file= file.path(output_dir, "cDC1&cDC2_ferrogene_information.csv"), row.names = FALSE)



### Find KEGG differential pathways
# Perform KEGG enrichment analysis on genes specific to cDC1 cells
entrezIDs <- bitr(cDC1_specific, fromType = "SYMBOL",
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = org.Hs.eg.db)
cDC1_gene <- entrezIDs$ENTREZID

cDC1_KEGG <- enrichKEGG(gene = cDC1_gene,
                        organism = "hsa",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)
cDC1_kk <- setReadable(cDC1_KEGG, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(cDC1_kk, file= file.path(output_dir, "cDC1&cDC2_cDC1_KEGG.csv"), quote = F)

pdf(file.path(output_dir, "cDC1&cDC2_cDC1_KEGG.pdf"), width = 9, height = 6)
dotplot(cDC1_KEGG, showCategory = 10, label_format = 100)
dev.off()


#Perform KEGG enrichment analysis on genes specific to cDC2 cells
entrezIDs <- bitr(cDC2_specific, fromType = "SYMBOL",
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = org.Hs.eg.db)
cDC2_gene <- entrezIDs$ENTREZID

cDC2_KEGG <- enrichKEGG(gene = cDC2_gene,
                        organism = "hsa",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)
cDC2_kk <- setReadable(cDC2_KEGG, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(cDC2_kk, file= file.path(output_dir, "cDC1&cDC2_cDC2_KEGG.csv"), quote = F)

pdf(file.path(output_dir, "cDC1&cDC2_cDC2_KEGG.pdf"), width = 9, height = 6)
dotplot(cDC2_KEGG, showCategory = 10, label_format = 100)
dev.off()




## Explain the differences in ferroptosis scores within the same cell subpopulation: CD8Trm vs. CD8Trm_exhausted
## Identify DEGs
# Load and Filter Differential Gene Data
cd8trm_diff <- read_excel("my_data/diff_genes_PAAD_immune_detailed_CD8Trm.xlsx")
cd8trm_exhausted_diff <- read_excel("my_data/diff_genes_PAAD_immune_detailed_CD8Trm_exhausted.xlsx")

cd8trm_diff$cell_type <- "CD8Trm"
cd8trm_exhausted_diff$cell_type <- "CD8Trm_exhausted"

combined_diff <- bind_rows(cd8trm_diff, cd8trm_exhausted_diff)
sig_genes <- combined_diff %>% 
  filter(PValue < 0.05 & abs(logFC) > 1)

#Obtain the specific differential genes of the two cell types
trm_specific <- setdiff(
  sig_genes$Gene[sig_genes$cell_type == "CD8Trm"],
  sig_genes$Gene[sig_genes$cell_type == "CD8Trm_exhausted"]
)

exhausted_specific <- setdiff(
  sig_genes$Gene[sig_genes$cell_type == "CD8Trm_exhausted"],
  sig_genes$Gene[sig_genes$cell_type == "CD8Trm"]
)

pdf(file.path(output_dir, "CD8Trm&CD8Trm_exhausted_diffgene_venn.pdf"), width = 14, height = 12)
venn.plot <- venn.diagram(
  x = list(
    CD8Trm = sig_genes$Gene[sig_genes$cell_type == "CD8Trm"],
    CD8Trm_exhausted = sig_genes$Gene[sig_genes$cell_type == "CD8Trm_exhausted"]
  ),
  filename = NULL,
  fill = c("#1b9e77", "#d95f02"),
  main = "Differential Genes Venn Diagram",
  main.cex = 2,
  cat.cex = 1.5,
  cex = 1.5
)
grid.draw(venn.plot)
dev.off()


#Output the specific differential genes of the two cell types and the information of common differential genes
common_genes <- intersect(
  sig_genes$Gene[sig_genes$cell_type == "CD8Trm"],
  sig_genes$Gene[sig_genes$cell_type == "CD8Trm_exhausted"]
)

create_sorted_df <- function(genes, cell_type) {
  df <- sig_genes %>% 
    filter(Gene %in% genes, cell_type == !!cell_type) %>%
    arrange(desc(abs(logFC))) %>%  
    select(Gene, logFC, PValue) %>%
    mutate(importance_rank = row_number())  
  return(df)
}

trm_specific_df <- create_sorted_df(trm_specific, "CD8Trm")
exhausted_specific_df <- create_sorted_df(exhausted_specific, "CD8Trm_exhausted")

common_df_trm <- create_sorted_df(common_genes, "CD8Trm") %>%
  rename(logFC_CD8Trm = logFC, PValue_CD8Trm = PValue)
common_df_exhausted <- create_sorted_df(common_genes, "CD8Trm_exhausted") %>%
  rename(logFC_CD8Trm_exhausted = logFC, PValue_CD8Trm_exhausted = PValue)

common_df <- common_df_trm %>%
  inner_join(common_df_exhausted, by = "Gene") %>%
  mutate(logFC_diff = abs(logFC_CD8Trm) - abs(logFC_CD8Trm_exhausted)) %>%
  arrange(desc(abs(logFC_CD8Trm) + abs(logFC_CD8Trm_exhausted)))  

wb <- createWorkbook()
addWorksheet(wb, "CD8Trm_specific")
addWorksheet(wb, "CD8Trm_exhausted_specific")
addWorksheet(wb, "Common_genes")
writeData(wb, "CD8Trm_specific", trm_specific_df)
writeData(wb, "CD8Trm_exhausted_specific", exhausted_specific_df)
writeData(wb, "Common_genes", common_df)

saveWorkbook(wb, file = file.path(output_dir, "CD8Trm&CD8Trm_exhausted_diffgene_sc.xlsx"), overwrite = TRUE)


#Obtain ferroptosis-related genes among the specific differential genes of the two cell types
trm_ferro_genes <- intersect(trm_specific, ferro_genes)
trm_ferro_ratio <- length(trm_ferro_genes) / length(trm_specific)

exhausted_ferro_genes <- intersect(exhausted_specific, ferro_genes)
exhausted_ferro_ratio <- length(exhausted_ferro_genes) / length(exhausted_specific)

ferro_results <- data.frame(
  CellType = c("CD8Trm", "CD8Trm_exhausted"),
  Total_Specific_Genes = c(length(trm_specific), length(exhausted_specific)),
  Ferroptosis_Genes = c(length(trm_ferro_genes), length(exhausted_ferro_genes)),
  Ferroptosis_Ratio = c(trm_ferro_ratio, exhausted_ferro_ratio),
  Ferroptosis_Gene_List = c(paste(trm_ferro_genes, collapse = "; "), 
                            paste(exhausted_ferro_genes, collapse = "; "))
)
write.csv(ferro_results, file= file.path(output_dir, "CD8Trm&CD8Trm_exhausted_ferrogene_information.csv"), row.names = FALSE)



### Find KEGG differential pathways
# Perform KEGG enrichment analysis on genes specific to CD8Trm cells
entrezIDs <- bitr(trm_specific, fromType = "SYMBOL",
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = org.Hs.eg.db)
cd8trm_gene <- entrezIDs$ENTREZID

cd8trm_KEGG <- enrichKEGG(gene = cd8trm_gene,
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
cd8trm_kk <- setReadable(cd8trm_KEGG, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(cd8trm_kk, file= file.path(output_dir, "CD8Trm&CD8Trm_exhausted_CD8Trm_KEGG.csv"), quote = F)

pdf(file.path(output_dir, "CD8Trm&CD8Trm_exhausted_CD8Trm_KEGG.pdf"), width = 9, height = 6)
dotplot(cd8trm_KEGG, showCategory = 10, label_format = 100)
dev.off()


#Perform KEGG enrichment analysis on genes specific to CD8Trm_exhausted cells
entrezIDs <- bitr(exhausted_specific, fromType = "SYMBOL",
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = org.Hs.eg.db)
cd8trm_exhausted_gene <- entrezIDs$ENTREZID

cd8trm_exhausted_KEGG <- enrichKEGG(gene = cd8trm_exhausted_gene,
                                    organism = "hsa",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.2,
                                    qvalueCutoff = 0.2)
cd8trm_exhausted_kk <- setReadable(cd8trm_exhausted_KEGG, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(cd8trm_exhausted_kk, file= file.path(output_dir, "CD8Trm&CD8Trm_exhausted_CD8Trmexhausted_KEGG.csv"), quote = F)

pdf(file.path(output_dir, "CD8Trm&CD8Trm_exhausted_CD8Trmexhausted_KEGG.pdf"), width = 9, height = 6)
dotplot(cd8trm_exhausted_KEGG, showCategory = 10, label_format = 100)
dev.off()




#### Mechanism of ferroptosis-sensitive cells nCAF
# Obtain GSEA significant pathways related to ferroptosis
diff_genes <- read_excel("my_data/diff_genes_PAAD_nonimmune_Fibroblast_nCAF.xlsx")

gene_rank <- diff_genes$logFC
names(gene_rank) <- diff_genes$Gene
gene_rank <- sort(gene_rank, decreasing = TRUE)

ferro_sets <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  filter(grepl("FERR|IRON|GLUTATHIONE|LIPID|FATTY_ACID|PEROXID|OXIDATIVE|GSH", 
               gs_name, ignore.case = TRUE))

set.seed(123) 
gsea_ferro <- GSEA(gene_rank, 
                   TERM2GENE = ferro_sets[, c("gs_name", "gene_symbol")],
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff = 1,
                   eps = 0,
                   seed = TRUE)
gsea_ferro_filtered <- gsea_ferro
gsea_ferro_filtered@result <- subset(gsea_ferro@result, pvalue < 0.05)

write.csv(as.data.frame(gsea_ferro_filtered), file = file.path(output_dir, "nCAF_GSEA_Ferroptosisresults.csv"), row.names = FALSE)



#Identify key driver genes
#Extract core genes from the significant pathways in the previous GSEA step
gsea_results <- read.csv("my_data/nCAF_GSEA_results.csv")
sig_pathways <- subset(gsea_results, pvalue < 0.05 )

core_genes <- unique(unlist(lapply(strsplit(sig_pathways$core_enrichment, "/"), function(x) {
  x[1:round(length(x)*0.2)]                          # Take the top 20% of genes in each pathway
})))

#Ferroptosis-related genes in the known FerrDb database
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
vertex_data_norm$pca_score <- pca_res$x[,1]       #Take the first principal component

#Save the list of key driver genes
hub_genes_pca <- vertex_data_norm %>%
  arrange(desc(pca_score))

write.csv(hub_genes_pca, file= file.path(output_dir, "nCAF_Ferroptosis_DriverGenes_Detailed.csv"), row.names = FALSE)

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

write.csv(driver_genes, file= file.path(output_dir, "nCAF_Ferroptosis_DriverGenes_Type.csv"),row.names = FALSE)

##Visual Network
pdf(file.path(output_dir, "nCAF_DriverGenes_PPI_Network.pdf"),  width = 14, height = 12)
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
  labs(title = "nCAF Ferroptosis Driver Genes PPI Network",
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

pdf(file.path(output_dir, "nCAF_DriverGenes_RadarPlot.pdf"), width = 10, height = 10)
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
#Read previously identified key driver genes
top_drivers <- hub_genes

#Read TRRUST database information
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

write.csv(tf_data_norm, file= file.path(output_dir, "nCAF_Ferroptosis_KeyTFs.csv"), row.names = FALSE)

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

pdf(file.path(output_dir, "nCAF_Ferroptosis_TF_Network.pdf"), width = 14, height = 12)
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
  labs(title = "nCAF Ferroptosis: Key Transcription Factors and Driver Genes") +
  theme_void() +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5))

dev.off()




### EMTEn vs. EMTEp
## Identify DEGs
# Filter Differential Gene Data
EMTEn_diff <- read_excel("my_data/diff_genes_PAAD_nonimmune_Endothelial_EMTEn.xlsx")
EMTEp_diff <- read_excel("my_data/diff_genes_PAAD_nonimmune_Epithelial_EMTEp.xlsx")

EMTEn_diff$cell_type <- "EMTEn"
EMTEp_diff$cell_type <- "EMTEp"

combined_diff <- bind_rows(EMTEn_diff, EMTEp_diff)

sig_genes <- combined_diff %>% 
  filter(PValue < 0.05 & abs(logFC) > 1)

#Obtain the specific differential genes of the two cell types
EMTEn_specific <- setdiff(
  sig_genes$Gene[sig_genes$cell_type == "EMTEn"],
  sig_genes$Gene[sig_genes$cell_type == "EMTEp"]
)

EMTEp_specific <- setdiff(
  sig_genes$Gene[sig_genes$cell_type == "EMTEp"],
  sig_genes$Gene[sig_genes$cell_type == "EMTEn"]
)

pdf(file.path(output_dir, "EMTEn&EMTEp_diffgene_venn.pdf"), width = 14, height = 12)
venn.plot <- venn.diagram(
  x = list(
    EMTEn = sig_genes$Gene[sig_genes$cell_type == "EMTEn"],
    EMTEp = sig_genes$Gene[sig_genes$cell_type == "EMTEp"]
  ),
  filename = NULL,
  fill = c("#1b9e77", "#d95f02"),
  main = "Differential Genes Venn Diagram",
  main.cex = 2,
  cat.cex = 1.5,
  cex = 1.5
)
grid.draw(venn.plot)
dev.off()

#Output the specific differential genes of the two cell types and the information of common differential genes
common_genes <- intersect(
  sig_genes$Gene[sig_genes$cell_type == "EMTEn"],
  sig_genes$Gene[sig_genes$cell_type == "EMTEp"]
)

create_sorted_df <- function(genes, cell_type) {
  df <- sig_genes %>% 
    filter(Gene %in% genes, cell_type == !!cell_type) %>%
    arrange(desc(abs(logFC))) %>%  
    select(Gene, logFC, PValue) %>%
    mutate(importance_rank = row_number())  
  
  return(df)
}

EMTEn_specific_df <- create_sorted_df(EMTEn_specific, "EMTEn")
EMTEp_specific_df <- create_sorted_df(EMTEp_specific, "EMTEp")

common_df_EMTEn <- create_sorted_df(common_genes, "EMTEn") %>%
  rename(logFC_EMTEn = logFC, PValue_EMTEn = PValue)
common_df_EMTEp <- create_sorted_df(common_genes, "EMTEp") %>%
  rename(logFC_EMTEp = logFC, PValue_EMTEp = PValue)

common_df <- common_df_EMTEn %>%
  inner_join(common_df_EMTEp, by = "Gene") %>%
  mutate(logFC_diff = abs(logFC_EMTEn) - abs(logFC_EMTEp)) %>%
  arrange(desc(abs(logFC_EMTEn) + abs(logFC_EMTEp)))  

wb <- createWorkbook()
addWorksheet(wb, "EMTEn_specific")
addWorksheet(wb, "EMTEp_specific")
addWorksheet(wb, "Common_genes")
writeData(wb, "EMTEn_specific", EMTEn_specific_df)
writeData(wb, "EMTEp_specific", EMTEp_specific_df)
writeData(wb, "Common_genes", common_df)

saveWorkbook(wb, file = file.path(output_dir, "EMTEn&EMTEp_diffgene_sc.xlsx"), overwrite = TRUE)

#Obtain ferroptosis-related genes among the specific differential genes of the two cell types
ferro_genes <- unique(ferr_db_genes$symbol)

EMTEn_ferro_genes <- intersect(EMTEn_specific, ferro_genes)
EMTEn_ferro_ratio <- length(EMTEn_ferro_genes) / length(EMTEn_specific)

EMTEp_ferro_genes <- intersect(EMTEp_specific, ferro_genes)
EMTEp_ferro_ratio <- length(EMTEp_ferro_genes) / length(EMTEp_specific)

ferro_results <- data.frame(
  CellType = c("EMTEn", "EMTEp"),
  Total_Specific_Genes = c(length(EMTEn_specific), length(EMTEp_specific)),
  Ferroptosis_Genes = c(length(EMTEn_ferro_genes), length(EMTEp_ferro_genes)),
  Ferroptosis_Ratio = c(EMTEn_ferro_ratio, EMTEp_ferro_ratio),
  Ferroptosis_Gene_List = c(paste(EMTEn_ferro_genes, collapse = "; "), 
                            paste(EMTEp_ferro_genes, collapse = "; "))
)
write.csv(ferro_results, file= file.path(output_dir, "EMTEn&EMTEp_ferrogene_information.csv"), row.names = FALSE)



### Find KEGG differential pathways
# Perform KEGG enrichment analysis on genes specific to EMTEn cells
entrezIDs <- bitr(EMTEn_specific, fromType = "SYMBOL",
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = org.Hs.eg.db)
EMTEn_gene <- entrezIDs$ENTREZID

EMTEn_KEGG <- enrichKEGG(gene = EMTEn_gene,
                         organism = "hsa",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)
EMTEn_kk <- setReadable(EMTEn_KEGG, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(EMTEn_kk, file = file.path(output_dir, "EMTEn&EMTEp_EMTEn_KEGG.csv"), quote = F)

pdf(file.path(output_dir, "EMTEn&EMTEn_EMTEn_KEGG.pdf"), width = 9, height = 6)
dotplot(EMTEn_KEGG, showCategory = 10, label_format = 100)
dev.off()

# Perform KEGG enrichment analysis on genes specific to EMTEp cells
entrezIDs <- bitr(EMTEp_specific, fromType = "SYMBOL",
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = org.Hs.eg.db)
EMTEp_gene <- entrezIDs$ENTREZID

EMTEp_KEGG <- enrichKEGG(gene = EMTEp_gene,
                         organism = "hsa",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)
EMTEp_kk <- setReadable(EMTEp_KEGG, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(EMTEp_kk, file = file.path(output_dir, "EMTEn&EMTEp_EMTEp_KEGG.csv"), quote = F)

pdf(file.path(output_dir, "EMTEn&EMTEp_EMTEp_KEGG.pdf"), width = 9, height = 6)
dotplot(EMTEp_KEGG, showCategory = 10, label_format = 100)
dev.off()






