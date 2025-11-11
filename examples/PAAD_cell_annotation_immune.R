######Immune cell annotation
######################################################



if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install(c("celldex","Seurat","dplyr","ggplot2","Matrix","tibble"))
library(Seurat)
library(dplyr)
library(ggplot2)
library(celldex)
library(Matrix)
library(tibble)



output_dir <- "my_results/immunecellannotation"



########################################### 1. Major Cell Population Identification
### 1.1 Cell Clustering
sce_umap <- readRDS("my_data/PAAD_immunecell_datapre.rds")
a<-as.matrix(GetAssayData(object = sce_umap@assays$RNA,layer="counts")[1:20,1:20])  
b<-as.matrix(GetAssayData(object = sce_umap@assays$RNA,layer="data")[1:20,1:20])    
c<-sce_umap@meta.data 

sce_cluster <- FindClusters(
  sce_umap,
  resolution = 2,  
  algorithm = 4,
  random.seed = 42
)
table(sce_cluster@meta.data$seurat_clusters)

DimPlot(sce_cluster, reduction = "umap", label = TRUE) + NoLegend()




### 1.2 Cell Annotation
# Find marker genes for all clusters
markers <- FindAllMarkers(
  sce_cluster,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 1,
  test.use = "wilcox"  
)

top10_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(markers,file= file.path(output_dir, "cluster_markers_PAAD_immune.csv"))
write.csv(top10_markers, file= file.path(output_dir, "Top10cluster_markers_PAAD_immune.csv"))

# Create marker gene sets
immune_markers <- list(
  Macrophage=c("CD14","CD68","CD163","LYZ","C1QC","CD4","KLF2"),
  DC=c("CLEC9A","CLEC10A","FLT3","SELL","CD4","CLEC4C"),
  Monocyte=c("CD14","LYZ","FCN1","S100A8","FCGR3A"),
  Mast=c("KIT","CPA3","TPSAB1","TPSB2","MS4A2"),
  NK=c("KLRD1","KLF2","TRDC","FCGR3A","NCAM1","NKG7","KLRF1","FGFBP2","CX3CR1"),
  T_cell=c("CD3D","CD3E","CD8A","CD2","CD4","CD7","CD52","CD3G"),
  B_cell=c("CD19","CD79A","MS4A1","CD79B","VPREB3")
)

FeaturePlot(sce_cluster, features = unlist(immune_markers), ncol = 4)
unique_markers <- unique(unlist(immune_markers))
DotPlot(sce_cluster, features = unique_markers) + 
  RotatedAxis() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red")
VlnPlot(sce_cluster, features = unique_markers, stack = T, flip = T) + 
  NoLegend()

# Annotation
new_cluster_ids <- c(
  "1" = "T",
  "2" = "T",
  "3" = "T",
  "4" = "Macrophage",
  "5" = "T",
  "6" = "Macrophage",
  "7" = "B",
  "8" = "T",
  "9" = "T",
  "10" = "DC",
  "11" = "Macrophage",
  "12" = "Monocyte",
  "13" = "NK",              
  "14" = "NK",  
  "15" = "Macrophage",     
  "16" = "NK",
  "17" = "Monocyte",
  "18" = "T",      
  "19" = "T",
  "20" = "T",
  "21" = "Mast",
  "22" = "B",
  "23" = "DC",
  "24" = "T"
)

sce_cluster$manual_annotation <- plyr::mapvalues(
  x = Idents(sce_cluster),
  from = names(new_cluster_ids),
  to = new_cluster_ids
)

pdf(file.path(output_dir, "UMAP_annotation_PAAD_immune.pdf"), width = 10, height = 8)
DimPlot(sce_cluster, group.by = "manual_annotation",label = TRUE,repel = TRUE) 
dev.off()

saveRDS(sce_cluster, "my_data/PAAD_immunecell_annotated_main.rds")




####################### 1.3 Convert barcode × gene matrix to celltype × gene matrix, aggregate expression values by cell type
expression_matrix <- GetAssayData(sce_cluster, assay = "RNA", slot = "counts")
expr_df <- as.data.frame(as.matrix(t(expression_matrix)))    # 转为barcode x gene
expr_df$celltype <- sce_cluster$manual_annotation

celltype_sum_matrix <- expr_df %>%
  group_by(celltype) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("celltype") %>%
  t()  

write.csv(celltype_sum_matrix, file= file.path(output_dir, "celltype_sum_countexpression_PAAD_immune.csv"))




################################################## 2. Further Subpopulation Analysis
### 2.1 T Cell Subpopulation Analysis
t_cells <- subset(sce_cluster, subset = manual_annotation == "T")

# Re-preprocessing
t_cells <- NormalizeData(t_cells, ormalization.method = "LogNormalize", scale.factor = 10000)
t_cells <- FindVariableFeatures(t_cells, selection.method = "vst", nfeatures = 3000)
t_cells <- ScaleData(t_cells)
t_cells <- RunPCA(t_cells, features = VariableFeatures(t_cells))
ElbowPlot(t_cells, ndims = 50)
cumulative_var <- cumsum(t_cells@reductions$pca@stdev^2 / sum(t_cells@reductions$pca@stdev^2))
elbow_point <- which.max(diff(cumulative_var) < 0.01)  
print(paste("Estimated elbow point at PC:", elbow_point))
plot(cumulative_var, type = "b", xlab = "PCs", ylab = "Cumulative Variance Explained")
abline(v = 25, col = "blue", lty = 2)
abline(h = 0.8, col = "grey", lty = 3)  
t_cells <- FindNeighbors(t_cells, reduction = "pca", dims = 1:25)
t_cells <- RunUMAP(t_cells, reduction = "pca", dims = 1:25)

# Clustering
t_cells <- FindClusters(t_cells, resolution = 2, algorithm = 4, random.seed = 42)
table(t_cells@meta.data$seurat_clusters)

DimPlot(t_cells, reduction = "umap", label = TRUE) + NoLegend()

# Find marker genes for T cell subpopulations
t_markers <- FindAllMarkers(
  t_cells,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 1,
  test.use = "wilcox"
)
top10_tmarkers <- t_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(t_markers, file= file.path(output_dir, "Tcell_cluster_markers_PAAD_immune.csv"))
write.csv(top10_tmarkers, file= file.path(output_dir, "Top10Tcell_cluster_markers_PAAD_immune.csv"))

# Create T cell marker sets
t_subset_markers <- list(
  CD8_Naive = c("CD8A", "CCR7", "LEF1", "TCF7", "SELL", "CD3D"),
  CD8_Tem_Effector = c("CD8A", "GZMB", "PRF1", "IFNG", "GNLY","GZMK", "CXCR3", "PTPRC", "CD3D"),
  CD8_Tcm = c("CD8A", "CCR7", "IL7R", "PTPRC", "CD3D"),   
  CD8_Trm = c("CD8A", "CD69", "ITGAE", "CXCR6", "ZNF683", "CD3D"),
  CD4_Naive = c("CD4", "CCR7", "LEF1", "TCF7", "SELL", "CD3D"),
  CD4_Tem_Effector = c("CD4", "GZMA", "IFNG", "CD3D","TNF", "GZMK", "CXCR3","PTPRC"),
  CD4_Tcm = c("CD4", "CCR7", "IL7R", "PTPRC", "CD3D"),
  CD4_Treg = c("CD4", "FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "CD3D"),
  GammaDelta = c("TRDC", "TRGV9", "TRDV2", "GZMB", "CD3D"),
  DN =c("TRAC","TRBC1","TRBC2","KLRB1","ZBTB16")
)

FeaturePlot(t_cells, features = unlist(t_subset_markers), ncol = 4)
DotPlot(t_cells, features = unique(unlist(t_subset_markers))) + 
  RotatedAxis() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red")
VlnPlot(t_cells, features = unique(unlist(t_subset_markers)), stack = T, flip = T) + 
  NoLegend()

# Annotation
new_t_labels <- c(
  "1" = "CD8 Tem/Teff",
  "2" = "CD4 Treg",
  "3" = "CD8 Trm",
  "4" = "CD8 Tem/Teff",
  "5" = "DN AlphaBeta T",
  "6" = "CD8 Tem/Teff",
  "7" = "CD8 Tcm",
  "8" = "Naive T",
  "9" = "CD4 Tcm",
  "10" = "CD8 Trm",
  "11" = "CD8 Trm_exhausted",
  "12" = "GammaDelta T",
  "13" = "Naive T",
  "14" = "DN AlphaBeta T",
  "15" = "CD8 Trm",
  "16" = "CD4 Tem/Teff",
  "17" = "CD8 Trm"
)

t_cells$t_subset_annotation <- plyr::mapvalues(
  x = Idents(t_cells),
  from = names(new_t_labels),
  to = new_t_labels
)

pdf(file.path(output_dir, "Tcell_UMAP_annotation_PAAD_immune.pdf"), width = 8, height = 6)
DimPlot(t_cells, group.by = "t_subset_annotation", label = TRUE, repel = TRUE)
dev.off()

saveRDS(t_cells, "my_data/PAAD_immunecell_annotated_Tcellsubset.rds")

# Convert T cell barcode × gene matrix to celltype × gene matrix for subsequent ferroptosis scoring
expression_matrix_Tcell <- GetAssayData(t_cells, assay = "RNA", slot = "counts")
expr_df_Tcell <- as.data.frame(as.matrix(t(expression_matrix_Tcell)))  
expr_df_Tcell$celltype <- t_cells$t_subset_annotation

celltype_sum_matrix_Tcell <- expr_df_Tcell %>%
  group_by(celltype) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("celltype") %>%
  t()  

write.csv(celltype_sum_matrix_Tcell, file= file.path(output_dir, "celltype_sum_countexpression_Tcellsubset_PAAD_immune.csv"))




###2.2 B Cell Subpopulation Analysis
B_cells <- subset(sce_cluster, subset = manual_annotation == "B")

#Re-preprocessing
B_cells <- NormalizeData(B_cells, ormalization.method = "LogNormalize", scale.factor = 10000)
B_cells <- FindVariableFeatures(B_cells, selection.method = "vst", nfeatures = 3000)
B_cells <- ScaleData(B_cells)
B_cells <- RunPCA(B_cells, features = VariableFeatures(B_cells))
ElbowPlot(B_cells, ndims = 50)
cumulative_var <- cumsum(B_cells@reductions$pca@stdev^2 / sum(B_cells@reductions$pca@stdev^2))
elbow_point <- which.max(diff(cumulative_var) < 0.01)  
print(paste("Estimated elbow point at PC:", elbow_point))
plot(cumulative_var, type = "b", xlab = "PCs", ylab = "Cumulative Variance Explained")
abline(v = 37, col = "blue", lty = 2)
abline(h = 0.8, col = "grey", lty = 3)  
B_cells <- FindNeighbors(B_cells, reduction = "pca", dims = 1:37)
B_cells <- RunUMAP(B_cells, reduction = "pca", dims = 1:37)

#Clustering
B_cells <- FindClusters(B_cells, resolution = 1, algorithm = 4, random.seed = 42)
table(B_cells@meta.data$seurat_clusters)

DimPlot(B_cells, reduction = "umap", label = TRUE) + NoLegend()

#Find marker genes for B cell subpopulations
B_markers <- FindAllMarkers(
  B_cells,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 1,
  test.use = "wilcox"
)
top10_Bmarkers <- B_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(B_markers, file= file.path(output_dir, "Bcell_cluster_markers_PAAD_immune.csv"))
write.csv(top10_Bmarkers, file= file.path(output_dir, "Top10Bcell_cluster_markers_PAAD_immune.csv"))

#Create B cell marker sets
B_subset_markers <- list(
  B_Naive = c("CD19","CD20","IGHD","IGHM","TCL1A","IL4R","CD27-", "CD38-"),
  B_Memory = c("CD27","CD21","TNFRSF13B","FCRL4","IGHG1", "IGHA1"),
  Plasma = c("CD38","CD138","PRDM1","XBP1","MZB1","IGHG3","JCHAIN")
)

FeaturePlot(B_cells, features = unlist(B_subset_markers), ncol = 4)
DotPlot(B_cells, features = unique(unlist(B_subset_markers))) + 
  RotatedAxis() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red")
VlnPlot(B_cells, features = unique(unlist(B_subset_markers)), stack = T, flip = T) + 
  NoLegend()

#Annotation
new_B_labels <- c(
  "1" = "Memory B",
  "2" = "Memory B",
  "3" = "Naive B",
  "4" = "Plasma"
)

B_cells$B_subset_annotation <- plyr::mapvalues(
  x = Idents(B_cells),
  from = names(new_B_labels),
  to = new_B_labels
)

pdf(file.path(output_dir, "Bcell_UMAP_annotation_PAAD_immune.pdf"), width = 8, height = 6)
DimPlot(B_cells, group.by = "B_subset_annotation", label = TRUE, repel = TRUE) 
dev.off()

saveRDS(B_cells, "my_data/PAAD_immunecell_annotated_Bcellsubset.rds")

#Convert B cell barcode × gene matrix to celltype × gene matrix for subsequent ferroptosis scoring
expression_matrix_Bcell <- GetAssayData(B_cells, assay = "RNA", slot = "counts")
expr_df_Bcell <- as.data.frame(as.matrix(t(expression_matrix_Bcell)))  
expr_df_Bcell$celltype <- B_cells$B_subset_annotation

celltype_sum_matrix_Bcell <- expr_df_Bcell %>%
  group_by(celltype) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("celltype") %>%
  t()  

write.csv(celltype_sum_matrix_Bcell, file= file.path(output_dir, "celltype_sum_countexpression_Bcellsubset_PAAD_immune.csv"))





###2.3 DC Cell Subpopulation Analysis
DC_cells <- subset(sce_cluster, subset = manual_annotation == "DC")

#Re-preprocessing
DC_cells <- NormalizeData(DC_cells, ormalization.method = "LogNormalize", scale.factor = 10000)
DC_cells <- FindVariableFeatures(DC_cells, selection.method = "vst", nfeatures = 3000)
DC_cells <- ScaleData(DC_cells)
DC_cells <- RunPCA(DC_cells, features = VariableFeatures(DC_cells))
ElbowPlot(DC_cells, ndims = 50)
cumulative_var <- cumsum(DC_cells@reductions$pca@stdev^2 / sum(DC_cells@reductions$pca@stdev^2))
elbow_point <- which.max(diff(cumulative_var) < 0.01)  
print(paste("Estimated elbow point at PC:", elbow_point))
plot(cumulative_var, type = "b", xlab = "PCs", ylab = "Cumulative Variance Explained")
abline(v = 34, col = "blue", lty = 2)
abline(h = 0.8, col = "grey", lty = 3)  
DC_cells <- FindNeighbors(DC_cells, reduction = "pca", dims = 1:34)
DC_cells <- RunUMAP(DC_cells, reduction = "pca", dims = 1:34)

#Clustering
DC_cells <- FindClusters(DC_cells, resolution = 1.5, algorithm = 4, random.seed = 42)
table(DC_cells@meta.data$seurat_clusters)

DimPlot(DC_cells, reduction = "umap", label = TRUE) + NoLegend()

#Find marker genes for DC cell subpopulations
DC_markers <- FindAllMarkers(
  DC_cells,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 1,
  test.use = "wilcox"
)
top10_DCmarkers <- DC_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(DC_markers, file= file.path(output_dir, "DCcell_cluster_markers_PAAD_immune.csv"))
write.csv(top10_DCmarkers, file= file.path(output_dir, "Top10DCcell_cluster_markers_PAAD_immune.csv"))

#Create DC cell marker sets
DC_subset_markers <- list(
  cDC1 = c("CLEC9A","XCR1","CADM1","PTDSS1","SCARB1","IL6ST","CD40","TNFRSF10B","IDO1","CST7","CLIC2","NET1","ANXA6","BATF3","ID2"),
  cDC2 = c("CD1C","CD1E","CLEC10A","FCER1A","FCGR2B","ADAM8","AXL","ADAM28","LY86","TIMM13","ARAF","HMGA1","PFDN1"),
  pDC = c("LILRA4","SLC15A4","PLD4","CCDC50","IL3RA","LY9","SELL","GAS6","TCF4","IRF7")
)

FeaturePlot(DC_cells, features = unlist(DC_subset_markers), ncol = 4)
DotPlot(DC_cells, features = unique(unlist(DC_subset_markers))) + 
  RotatedAxis() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red")
VlnPlot(DC_cells, features = unique(unlist(DC_subset_markers)), stack = T, flip = T) + 
  NoLegend()

#Annotation
new_DC_labels <- c(
  "1" = "cDC2",
  "2" = "cDC2",
  "3" = "cDC2",
  "4" = "pDC",
  "5" = "cDC1",
  "6" = "cDC2"
)

DC_cells$DC_subset_annotation <- plyr::mapvalues(
  x = Idents(DC_cells),
  from = names(new_DC_labels),
  to = new_DC_labels
)

pdf(file.path(output_dir, "DCcell_UMAP_annotation_PAAD_immune.pdf"), width = 8, height = 6)
DimPlot(DC_cells, group.by = "DC_subset_annotation", label = TRUE, repel = TRUE) 
dev.off()

saveRDS(DC_cells, "my_data/PAAD_immunecell_annotated_DCcellsubset.rds")

#Convert DC cell barcode × gene matrix to celltype × gene matrix for subsequent ferroptosis scoring
expression_matrix_DCcell <- GetAssayData(DC_cells, assay = "RNA", slot = "counts")
expr_df_DCcell <- as.data.frame(as.matrix(t(expression_matrix_DCcell)))  
expr_df_DCcell$celltype <- DC_cells$DC_subset_annotation

celltype_sum_matrix_DCcell <- expr_df_DCcell %>%
  group_by(celltype) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("celltype") %>%
  t()  

write.csv(celltype_sum_matrix_DCcell, file= file.path(output_dir, "celltype_sum_countexpression_DCcellsubset_PAAD_immune.csv"))




############################################### 3. Detailed Cell Population Identification
# Merge T cell, B cell, and DC cell subpopulation information back to original data
all_barcodes <- colnames(sce_cluster)
detailed_annotations <- as.character(sce_cluster$manual_annotation)

t_cell_mask <- all_barcodes %in% colnames(t_cells)
b_cell_mask <- all_barcodes %in% colnames(B_cells)
dc_cell_mask <- all_barcodes %in% colnames(DC_cells)

t_labels <- as.character(t_cells$t_subset_annotation)[match(all_barcodes[t_cell_mask], colnames(t_cells))]
b_labels <- as.character(B_cells$B_subset_annotation)[match(all_barcodes[b_cell_mask], colnames(B_cells))]
dc_labels <- as.character(DC_cells$DC_subset_annotation)[match(all_barcodes[dc_cell_mask], colnames(DC_cells))]

detailed_annotations[t_cell_mask] <- t_labels
detailed_annotations[b_cell_mask] <- b_labels
detailed_annotations[dc_cell_mask] <- dc_labels

celltype_levels <- c("Macrophage","Monocyte","Mast","NK",
                     "cDC1","cDC2","pDC",
                     "Naive B","Memory B","Plasma",
                     "CD8 Tem/Teff","CD8 Tcm","CD8 Trm",
                     "CD4 Tem/Teff","CD4 Tcm","CD4 Treg",
                     "Naive T","GammaDelta T","DN AlphaBeta T",
                     "CD8 Trm_exhausted")

sce_cluster$detailed_annotation <- factor(detailed_annotations, levels = celltype_levels)
table(sce_cluster$detailed_annotation, useNA = "always")
saveRDS(sce_cluster, "my_data/PAAD_immunecell_annotated_detailed.rds")

#Results after merging annotations
celltype_colors <- c(
  "Macrophage" = "#1F77B4",
  "Monocyte" = "#AEC7E8",
  "Mast" = "#FF7F0E",
  "NK" = "#FFBB78",
  "cDC1" = "#2CA02C",
  "cDC2" = "#98DF8A",
  "pDC" = "#D62728",
  "Naive B" = "#9467BD",
  "Memory B" = "#C5B0D5",
  "Plasma" = "#8C564B",
  "CD8 Tem/Teff" = "#E0D1D7",
  "CD8 Tcm" = "#F7B6D2",
  "CD8 Trm" = "#E377C2",
  "CD4 Tem/Teff" = "#6B6ECF",
  "CD4 Tcm" = "#393B79",
  "CD4 Treg" = "#D0CCF0",
  "Naive T" = "#DBDB8D",
  "GammaDelta T" = "#17BECF",
  "DN AlphaBeta T" = "#9EDAE5",
  "CD8 Trm_exhausted" = "#7F7F7F"
)

pdf(file.path(output_dir, "UMAP_annotation_detailed_PAAD_immune.pdf"), width = 12, height = 8)
DimPlot(sce_cluster, 
        group.by = "detailed_annotation",
        label = TRUE,
        repel = TRUE,
        cols = celltype_colors) +
  ggtitle("PAAD Immune Cell Annotation") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

#Convert merged barcode × gene matrix to celltype × gene matrix for subsequent ferroptosis score calculation
expression_matrix_detailed <- GetAssayData(sce_cluster, assay = "RNA", slot = "counts")
expr_df_detailed <- as.data.frame(as.matrix(t(expression_matrix_detailed)))  
expr_df_detailed$celltype <- sce_cluster$detailed_annotation

celltype_sum_matrix_detailed <- expr_df_detailed %>%
  group_by(celltype) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("celltype") %>%
  t()  

write.csv(celltype_sum_matrix_detailed, file= file.path(output_dir, "celltype_sum_expression_detailed.csv"))



