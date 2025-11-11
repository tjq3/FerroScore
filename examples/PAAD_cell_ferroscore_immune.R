######Calculate the ferroptosis score of each type of immune cell
################################################################




library(FerroScore)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("edgeR", "stringr", "igraph", "ggraph", "STRINGdb",
                       "dplyr", "clusterProfiler", "org.Hs.eg.db","enrichplot",
                       "ggplot2", "GOplot", "openxlsx","RColorBrewer"))

library(edgeR)
library(stringr)
library(igraph)
library(ggraph)
library(STRINGdb)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(GOplot)
library(openxlsx)
library(RColorBrewer)


output_dir <- "my_results/immunecellferroscore"




#Calculate the ferroptosis score for each type of immune cell individually (taking Macro as an example)
counts_data <- read.csv("my_data/celltype_count_PAAD_immune_Macro.csv", row.names = 1)
meta_data <- read.csv("my_data/celltype_meta_PAAD_immune_Macro.csv", stringsAsFactors = TRUE)

result_Macro<- run_ferroptosis_analysis(
  counts_data = counts_data,
  meta_data = meta_data,
  control_group = "control",
  species = "human",
  ppi_file = NULL,
  output_dir = tempdir(),
  save_intermediate = FALSE,
  parallel = FALSE,
  n_cores = 1)

result_Macro$summary
result_Macro$merged_DEG
result_Macro$results

ferroptosis_scores <- t(result_Macro$summary)
write.csv(ferroptosis_scores, file= file.path(output_dir, "ferroptosis_scores_PAAD_immune_Macro.csv"), row.names = TRUE)


deg_list <- split(result_Macro$merged_DEG, result_Macro$merged_DEG$Sample)
names(deg_list)
wb <- createWorkbook()
for (sample_name in names(deg_list)) {
  addWorksheet(wb, sheetName = sample_name)
  writeData(wb, sheet = sample_name, x = deg_list[[sample_name]])
}
saveWorkbook(wb, file = file.path(output_dir, "diff_genes_PAAD_immune_Macro.xlsx"), overwrite = TRUE)




#Pie chart of immune cells, showing the number of each type of immune cell
sce_cluster <- readRDS("my_data/PAAD_immunecell_annotated_detailed.rds")
cell_counts <- table(sce_cluster$detailed_annotation)
count_df <- as.data.frame(cell_counts)
colnames(count_df) <- c("CellType", "Count")
ordered_colors <- celltype_colors[names(celltype_colors) %in% count_df$CellType]


pdf(file.path(output_dir, "Celltype_piechart.pdf"), width = 10, height = 8)
ggplot(count_df, aes(x = "", y = Count, fill = CellType)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = ordered_colors) +
  theme_void() +
  labs(title = "PAAD Immune Cell Composition") +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "right",
        legend.title = element_blank()) +
  geom_text(aes(label = paste0(CellType, "\n", Count)), 
            position = position_stack(vjust = 0.5),
            size = 3)
dev.off()




#Bar Chart of ferroptosis score
ferroptosis_data_PAADimmune <- data.frame(
  CellType = c("Macrophage", "Monocyte", "Mast", "NK", "cDC1", "cDC2", "pDC", 
               "Naive B", "Memory B", "Plasma", "CD8 Tem/Teff", "CD8 Tcm", 
               "CD8 Trm", "CD4 Tem/Teff", "CD4 Tcm", "CD4 Treg", "Naive T", 
               "GammaDelta T", "DN AlphaBeta T", "CD8 Trm_exhausted"),
  FerroptosisScore = c(1.436645, 0.4958272, 0.4800443, 0.7276803, 0.6876329, 
                       0.5618609, NA, NA, 0.6005962, NA, NA, 1.081676, 
                       0.7132325, NA, NA, NA, NA, 0.6686389, NA, 0.4824523)
)
ferroptosis_data_PAADimmune <- ferroptosis_data_PAADimmune[!is.na(ferroptosis_data_PAADimmune$FerroptosisScore), ]

ferroptosis_data_PAADimmune <- ferroptosis_data_PAADimmune[order(-ferroptosis_data_PAADimmune$FerroptosisScore), ]
ferroptosis_data_PAADimmune$CellType <- factor(ferroptosis_data_PAADimmune$CellType, 
                                               levels = ferroptosis_data_PAADimmune$CellType)


pdf(file.path(output_dir, "Ferroptosisscore_barplot.pdf"), width = 10, height = 6)
ggplot(ferroptosis_data_PAADimmune, aes(x = CellType, y = FerroptosisScore, fill = FerroptosisScore)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#D8BFD8", high = "#4B0082") + 
  labs(x = "Cell Type",
       y = "Ferroptosis Score") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.8),
    axis.ticks = element_line(colour = "black", linewidth = 0.8), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(face = "bold", size = 18),
    legend.position = "right",
    panel.background = element_blank()) +
  geom_text(aes(label = round(FerroptosisScore, 3)), vjust = -0.3, size = 3)
dev.off()




#UMAP of ferroptosis index
celltype_to_score <- c(
  "Macrophage" = 1.436645,
  "Monocyte" = 0.4958272,
  "Mast" = 0.4800443,
  "NK" = 0.7276803,
  "cDC1" = 0.6876329,
  "cDC2" = 0.5618609,
  "pDC" = NA,
  "Naive B" = NA,
  "Memory B" = 0.6005962,
  "Plasma" = NA,
  "CD8 Tem/Teff" = NA,
  "CD8 Tcm" = 1.081676,
  "CD8 Trm" = 0.7132325,
  "CD4 Tem/Teff" = NA,
  "CD4 Tcm" = NA,
  "CD4 Treg" = NA,
  "Naive T" = NA,
  "GammaDelta T" = 0.6686389,
  "DN AlphaBeta T" = NA,
  "CD8 Trm_exhausted" = 0.4824523
)

ferro_scores <- rep(NA, ncol(sce_cluster))
names(ferro_scores) <- colnames(sce_cluster)
cell_types <- as.character(sce_cluster$detailed_annotation)

for(ctype in names(celltype_to_score)){
  idx <- which(cell_types == ctype)
  if(length(idx) > 0){
    ferro_scores[idx] <- celltype_to_score[ctype]
  }
}

sce_cluster$ferroptosis_score <- ferro_scores[colnames(sce_cluster)]


calculate_ferroptosis_index <- function(scores) {
  valid_scores <- scores[!is.na(scores)]
  quantile_cutoffs <- quantile(valid_scores, probs = seq(0, 1, 0.1))
  assign_index <- function(x) {
    if(is.na(x)) return(1)  
    for(i in 10:1) {
      if(x >= quantile_cutoffs[i]) return(i)
    }
    return(1)
  }
  indices <- sapply(scores, assign_index)
  return(indices)
}

sce_cluster$ferroptosis_index <- calculate_ferroptosis_index(sce_cluster$ferroptosis_score)


ferroptosis_colors <- colorRampPalette(c("#F2E0F7", "#4B0082"))(10)  
names(ferroptosis_colors) <- 1:10
pdf(file.path(output_dir, "Ferroptosisindex_UMAP.pdf"), width = 12, height = 8)
DimPlot(sce_cluster, 
        group.by = "ferroptosis_index",
        cols = ferroptosis_colors,
        order = order(as.numeric(sce_cluster$ferroptosis_index)),  
        pt.size = 0.5) +
  ggtitle("PAAD Immune Cells by Ferroptosis Index") +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "right") +
  guides(color = guide_legend(title = "Ferroptosis\nIndex", 
                              override.aes = list(size = 4)))
dev.off()


print(table(sce_cluster$ferroptosis_index, useNA = "always"))
print(aggregate(ferroptosis_index ~ detailed_annotation, data = sce_cluster@meta.data, 
                FUN = function(x) round(mean(x, na.rm = TRUE), 2)))
saveRDS(sce_cluster, "my_data/PAAD_immunecell_with_ferroptosis.rds")






