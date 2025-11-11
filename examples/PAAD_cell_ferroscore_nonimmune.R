######Calculate the ferroptosis score of each type of nonimmune cell
###################################################################



library(FerroScore)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("edgeR", "stringr", "igraph", "ggraph", "STRINGdb",
                       "dplyr", "clusterProfiler", "org.Hs.eg.db",
                       "enrichplot", "ggplot2", "GOplot", "openxlsx","RColorBrewer"))

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



output_dir <- "my_results/nonimmunecellferroscore"



##Calculate the ferroptosis score for each type of nonimmune cell individually (taking Endothelial as an example)
counts_data <- read.csv("my_data/celltype_count_PAAD_nonimmune_Endothelial.csv", row.names = 1)
meta_data <- read.csv("my_data/celltype_meta_PAAD_nonimmune_Endothelial.csv", stringsAsFactors = TRUE)

result_Endothelial<- run_ferroptosis_analysis(
  counts_data = counts_data,
  meta_data = meta_data,
  control_group = "control",
  species = "human",
  ppi_file = NULL,
  output_dir = tempdir(),
  save_intermediate = FALSE,
  parallel = FALSE,
  n_cores = 2)

result_Endothelial$summary
result_Endothelial$merged_DEG
result_Endothelial$results

ferroptosis_scores <- t(result_Endothelial$summary)
write.csv(ferroptosis_scores, file = file.path(output_dir, "ferroptosis_scores_PAAD_nonimmune_Endothelial.csv"), row.names = TRUE)


deg_list <- split(result_Endothelial$merged_DEG, result_Endothelial$merged_DEG$Sample)
names(deg_list)
wb <- createWorkbook()
for (sample_name in names(deg_list)) {
  addWorksheet(wb, sheetName = sample_name)
  writeData(wb, sheet = sample_name, x = deg_list[[sample_name]])
}
saveWorkbook(wb, file = file.path(output_dir, "diff_genes_PAAD_nonimmune_Endothelial.xlsx"), overwrite = TRUE)




#Pie chart of five major nonimmune cell types
sce_cluster <- readRDS("my_data/PAAD_nonimmunecell_annotated_detailed.rds")
cell_counts <- table(sce_cluster$manual_annotation)
count_df <- as.data.frame(cell_counts)
colnames(count_df) <- c("CellType", "Count")
ordered_colors <- celltype_colors[names(celltype_colors) %in% count_df$CellType]

pdf(file.path(output_dir, "Celltype_piechart_main.pdf"), width = 10, height = 8)
ggplot(count_df, aes(x = "", y = Count, fill = CellType)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = ordered_colors) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "right",
        legend.title = element_blank()) +
  geom_text(aes(label = paste0(CellType, "\n", Count)), 
            position = position_stack(vjust = 0.5),
            size = 3)
dev.off()




#UMAP of ferroptosis index in five major nonimmune cell types
celltype_to_score <- c(
  "Epithelial" = 0.455144,
  "Endothelial" = 1.249963,
  "Fibroblast" = NA,
  "Stellate" = 0.557225,
  "Schwann" = NA
)

ferro_scores <- rep(NA, ncol(sce_cluster))
names(ferro_scores) <- colnames(sce_cluster)
cell_types <- as.character(sce_cluster$manual_annotation)

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
pdf(file.path(output_dir, "Ferroptosisindex_UMAP_main.pdf"), width = 12, height = 8)
DimPlot(sce_cluster, 
        group.by = "ferroptosis_index",
        cols = ferroptosis_colors,
        order = order(as.numeric(sce_cluster$ferroptosis_index)),  
        pt.size = 0.5) +
  ggtitle("PAAD Nonimmune Cells by Ferroptosis Index") +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "right") +
  guides(color = guide_legend(title = "Ferroptosis\nIndex", 
                              override.aes = list(size = 4)))
dev.off()


print(table(sce_cluster$ferroptosis_index, useNA = "always"))
print(aggregate(ferroptosis_index ~ manual_annotation, data = sce_cluster@meta.data, 
                FUN = function(x) round(mean(x, na.rm = TRUE), 2)))
saveRDS(sce_cluster, "my_data/PAAD_nonimmunecell_main_with_ferroptosis.rds")




#After subdividing Epithelial, Endothelial, Fibroblast, and Stellate subtypes, UMAP of ferroptosis index for each cell type
#Epithelial
Epithelial_cells <- readRDS("my_data/PAAD_nonimmunecell_annotated_Epithelialsubset.rds")
celltype_to_score_Ep <- c(
  "EMTEp" = 0.5819942,
  "metEp" = NA,
  "imEp" = 0.6421859,
  "pEp" = NA,
  "sEp" = NA,
  "ciliumEp" = 0.4848696,
  "sarcomatoidEp" = NA
)

ferro_scores_Ep <- rep(NA, ncol(Epithelial_cells))
names(ferro_scores_Ep) <- colnames(Epithelial_cells)
cell_types_Ep <- as.character(Epithelial_cells$Epithelial_subset_annotation)
for(ctype in names(celltype_to_score_Ep)){
  idx <- which(cell_types_Ep == ctype)
  if(length(idx) > 0){
    ferro_scores_Ep[idx] <- celltype_to_score_Ep[ctype]
  }
}
Epithelial_cells$ferroptosis_score <- ferro_scores_Ep[colnames(Epithelial_cells)]


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
Epithelial_cells$ferroptosis_index <- calculate_ferroptosis_index(Epithelial_cells$ferroptosis_score)


ferroptosis_colors <- colorRampPalette(c("#F2E0F7", "#4B0082"))(10)  
names(ferroptosis_colors) <- 1:10
pdf(file.path(output_dir, "Ferroptosisindex_UMAP_Epithelial.pdf"), width = 12, height = 8)
DimPlot(Epithelial_cells, 
        group.by = "ferroptosis_index",
        cols = ferroptosis_colors,
        order = order(as.numeric(Epithelial_cells$ferroptosis_index)),  
        pt.size = 0.5) +
  ggtitle("Epithelial Subset Cell by Ferroptosis Index") +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "right") +
  guides(color = guide_legend(title = "Ferroptosis\nIndex", 
                              override.aes = list(size = 4)))
dev.off()

print(table(Epithelial_cells$ferroptosis_index, useNA = "always"))
print(aggregate(ferroptosis_index ~ Epithelial_subset_annotation, data = Epithelial_cells@meta.data, 
                FUN = function(x) round(mean(x, na.rm = TRUE), 2)))
saveRDS(Epithelial_cells, "my_data/PAAD_nonimmunecell_Epithelial_with_ferroptosis.rds")




##Fibroblast
Fibroblast_cells <- readRDS("my_data/PAAD_nonimmunecell_annotated_Fibroblastsubset.rds")
celltype_to_score_CAF <- c(
  "mCAF" = 0.5344095,
  "iCAF" = 0.6401759,
  "hspCAF" = 0.6679029,
  "dCAF" = NA,
  "nCAF" = 2.073027
)

ferro_scores_CAF <- rep(NA, ncol(Fibroblast_cells))
names(ferro_scores_CAF) <- colnames(Fibroblast_cells)
cell_types_CAF <- as.character(Fibroblast_cells$Fibroblast_subset_annotation)
for(ctype in names(celltype_to_score_CAF)){
  idx <- which(cell_types_CAF == ctype)
  if(length(idx) > 0){
    ferro_scores_CAF[idx] <- celltype_to_score_CAF[ctype]
  }
}
Fibroblast_cells$ferroptosis_score <- ferro_scores_CAF[colnames(Fibroblast_cells)]


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
Fibroblast_cells$ferroptosis_index <- calculate_ferroptosis_index(Fibroblast_cells$ferroptosis_score)


ferroptosis_colors <- colorRampPalette(c("#F2E0F7", "#4B0082"))(10)  
names(ferroptosis_colors) <- 1:10
pdf(file.path(output_dir, "Ferroptosisindex_UMAP_Fibroblast.pdf"), width = 12, height = 8)
DimPlot(Fibroblast_cells, 
        group.by = "ferroptosis_index",
        cols = ferroptosis_colors,
        order = order(as.numeric(Fibroblast_cells$ferroptosis_index)),  
        pt.size = 0.5) +
  ggtitle("Fibroblast Subset Cell by Ferroptosis Index") +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "right") +
  guides(color = guide_legend(title = "Ferroptosis\nIndex", 
                              override.aes = list(size = 4)))
dev.off()

print(table(Fibroblast_cells$ferroptosis_index, useNA = "always"))
print(aggregate(ferroptosis_index ~ Fibroblast_subset_annotation, data = Fibroblast_cells@meta.data, 
                FUN = function(x) round(mean(x, na.rm = TRUE), 2)))
saveRDS(Fibroblast_cells, "my_data/PAAD_nonimmunecell_Fibroblast_with_ferroptosis.rds")




#Stellate
Stellate_cells <- readRDS("my_data/PAAD_nonimmunecell_annotated_Stellatesubset.rds")
celltype_to_score_PSC <- c(
  "fPSC" = NA,
  "metPSC"= NA,
  "iPSC" = 0.8332458,
  "imPSC" = NA,
  "nPSC" = 0.5341263,
  "cPSC" = 0.8013463
)

ferro_scores_PSC <- rep(NA, ncol(Stellate_cells))
names(ferro_scores_PSC) <- colnames(Stellate_cells)
cell_types_PSC <- as.character(Stellate_cells$Stellate_subset_annotation)
for(ctype in names(celltype_to_score_PSC)){
  idx <- which(cell_types_PSC == ctype)
  if(length(idx) > 0){
    ferro_scores_PSC[idx] <- celltype_to_score_PSC[ctype]
  }
}
Stellate_cells$ferroptosis_score <- ferro_scores_PSC[colnames(Stellate_cells)]


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
Stellate_cells$ferroptosis_index <- calculate_ferroptosis_index(Stellate_cells$ferroptosis_score)


ferroptosis_colors <- colorRampPalette(c("#F2E0F7", "#4B0082"))(10)  
names(ferroptosis_colors) <- 1:10
pdf(file.path(output_dir, "Ferroptosisindex_UMAP_Stellate.pdf"), width = 12, height = 8)
DimPlot(Stellate_cells, 
        group.by = "ferroptosis_index",
        cols = ferroptosis_colors,
        order = order(as.numeric(Stellate_cells$ferroptosis_index)),  
        pt.size = 0.5) +
  ggtitle("Stellate Subset Cell by Ferroptosis Index") +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "right") +
  guides(color = guide_legend(title = "Ferroptosis\nIndex", 
                              override.aes = list(size = 4)))
dev.off()

print(table(Stellate_cells$ferroptosis_index, useNA = "always"))
print(aggregate(ferroptosis_index ~ Stellate_subset_annotation, data = Stellate_cells@meta.data, 
                FUN = function(x) round(mean(x, na.rm = TRUE), 2)))
saveRDS(Stellate_cells, "my_data/PAAD_nonimmunecell_Stellate_with_ferroptosis.rds")




#Endothelial
Endothelial_cells <- readRDS("my_data/PAAD_nonimmunecell_annotated_Endothelialsubset.rds")
celltype_to_score_En <- c(
  "lEn" = NA,
  "vEn_Artery"= NA,
  "vEn_Vein" = NA,
  "vEn_Capillary" = NA,
  "iEn" = 0.4983876,
  "EMTEn" = 0.6732866
)

ferro_scores_En <- rep(NA, ncol(Endothelial_cells))
names(ferro_scores_En) <- colnames(Endothelial_cells)
cell_types_En <- as.character(Endothelial_cells$Endothelial_subset_annotation)
for(ctype in names(celltype_to_score_En)){
  idx <- which(cell_types_En == ctype)
  if(length(idx) > 0){
    ferro_scores_En[idx] <- celltype_to_score_En[ctype]
  }
}
Endothelial_cells$ferroptosis_score <- ferro_scores_En[colnames(Endothelial_cells)]


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
Endothelial_cells$ferroptosis_index <- calculate_ferroptosis_index(Endothelial_cells$ferroptosis_score)


ferroptosis_colors <- colorRampPalette(c("#F2E0F7", "#4B0082"))(10)  
names(ferroptosis_colors) <- 1:10
pdf(file.path(output_dir, "Ferroptosisindex_UMAP_Endothelial.pdf"), width = 12, height = 8)
DimPlot(Endothelial_cells, 
        group.by = "ferroptosis_index",
        cols = ferroptosis_colors,
        order = order(as.numeric(Endothelial_cells$ferroptosis_index)),  
        pt.size = 0.5) +
  ggtitle("Endothelial Subset Cell by Ferroptosis Index") +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "right") +
  guides(color = guide_legend(title = "Ferroptosis\nIndex", 
                              override.aes = list(size = 4)))
dev.off()

print(table(Endothelial_cells$ferroptosis_index, useNA = "always"))
print(aggregate(ferroptosis_index ~ Endothelial_subset_annotation, data = Endothelial_cells@meta.data, 
                FUN = function(x) round(mean(x, na.rm = TRUE), 2)))
saveRDS(Endothelial_cells, "my_data/PAAD_nonimmunecell_Endothelial_with_ferroptosis.rds")




#Total Bar Chart of ferroptosis score
ferroptosis_data_PAADnonimmune <- data.frame(
  CellType = c("EMTEp", "metEp", "imEp", "pEp", "sEp", "ciliumEp", "sarcomatoidEp", 
               "mCAF", "iCAF", "hspCAF", "dCAF", "nCAF", 
               "fPSC", "metPSC", "iPSC", "imPSC", "nPSC","cPSC", 
               "lEn", "vEn_Artery", "vEn_Vein","vEn_Capillary","iEn","EMTEn"),
  FerroptosisScore = c(0.5819942, NA, 0.6421859, NA, NA,0.4848696,NA, 
                       0.5344095, 0.6401759, 0.6679029, NA, 2.073027, 
                       NA, NA, 0.8332458,NA,0.5341263,0.8013463,
                       NA,NA,NA,NA,0.4983876,0.6732866)
)

ferroptosis_data_PAADnonimmune$Group <- NA
ferroptosis_data_PAADnonimmune$Group[grep("Ep$|Ep", ferroptosis_data_PAADnonimmune$CellType)] <- "Ep"
ferroptosis_data_PAADnonimmune$Group[grep("CAF$|CAF", ferroptosis_data_PAADnonimmune$CellType)] <- "CAF"
ferroptosis_data_PAADnonimmune$Group[grep("PSC$|PSC", ferroptosis_data_PAADnonimmune$CellType)] <- "PSC"
ferroptosis_data_PAADnonimmune$Group[grep("En$|En", ferroptosis_data_PAADnonimmune$CellType)] <- "En"

ferroptosis_data <- ferroptosis_data_PAADnonimmune[!is.na(ferroptosis_data_PAADnonimmune$FerroptosisScore), ]
ferroptosis_data <- ferroptosis_data[order(ferroptosis_data$Group, -ferroptosis_data$FerroptosisScore), ]
ferroptosis_data$CellType <- factor(ferroptosis_data$CellType, levels = ferroptosis_data$CellType)


pdf(file.path(output_dir, "Ferroptosisscore_barplotall.pdf"), width = 12, height = 6)
ggplot(ferroptosis_data, aes(x = CellType, y = FerroptosisScore, fill = FerroptosisScore)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#D8BFD8", high = "#4B0082", 
                      name = "Ferroptosis Score") + 
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



