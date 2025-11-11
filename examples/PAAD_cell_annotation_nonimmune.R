######Nonimmune cell annotation
######################################################




if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('celldex')
library(Seurat)
library(dplyr)
library(ggplot2)
library(celldex)
library(Matrix)
library(tibble)
library(clusterProfiler)


output_dir <- "my_results/nonimmunecellannotation"





########################################### 1. Major Cell Population Identification
### 1.1 Cell Clustering
sce_umap <- readRDS("my_data/PAAD_nonimmunecell_datapre.rds")
a <- as.matrix(sce_umap[["RNA"]]$counts[1:20, 1:20])
b <- as.matrix(sce_umap[["RNA"]]$data[1:20, 1:20])
c<-sce_umap@meta.data 

sce_cluster <- FindClusters(
  sce_umap,
  resolution = 1,  
  algorithm = 4,
  random.seed = 42
)
table(sce_cluster@meta.data$seurat_clusters)

DimPlot(sce_cluster, reduction = "umap", label = TRUE) + NoLegend()



### 1.2 Cell Annotation
# Find marker genes for all clusters
sce_cluster <- JoinLayers(sce_cluster)

markers <- FindAllMarkers(
  sce_cluster,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 1,
  test.use = "wilcox"  
)

top20_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

write.csv(markers, file= file.path(output_dir, "cluster_markers_PAAD_nonimmune.csv"))
write.csv(top20_markers, file= file.path(output_dir, "Top20cluster_markers_PAAD_nonimmune.csv"))

# Create marker gene sets
nonimmune_markers <- list(
  Epithelial = c("KRT8","KRT18","KRT19","CDH1","EPCAM","CLDN4"),
  Endothelial = c("PECAM1","CDH5","CD34","NOS3","VWF"),
  Fibroblast = c("FAP","PDGFRA","FGF7","LUM","DCN"),
  Stellate = "RGS5",
  Schwann = "SOX10"
)

FeaturePlot(sce_cluster, features = unlist(nonimmune_markers), ncol = 4)
unique_markers <- unique(unlist(nonimmune_markers))
DotPlot(sce_cluster, features = unique_markers) + 
  RotatedAxis() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red")
VlnPlot(sce_cluster, features = unique_markers, stack = T, flip = T) + 
  NoLegend()

# Annotation
new_cluster_ids <- c(
  "1" = "Fibroblast",
  "2" = "Stellate",
  "3" = "Fibroblast",
  "4" = "Epithelial",
  "5" = "Fibroblast",
  "6" = "Endothelial",
  "7" = "Fibroblast",
  "8" = "Epithelial",
  "9" = "Fibroblast",
  "10" = "Epithelial",
  "11" = "Fibroblast",
  "12" = "Epithelial",
  "13" = "Fibroblast",              
  "14" = "Epithelial",  
  "15" = "Epithelial",     
  "16" = "Epithelial",
  "17" = "Epithelial",
  "18" = "Epithelial",      
  "19" = "Endothelial",
  "20" = "Epithelial",
  "21" = "Epithelial",
  "22" = "Epithelial",
  "23" = "Epithelial",
  "24" = "Schwann",
  "25" = "Fibroblast"
)

sce_cluster$manual_annotation <- plyr::mapvalues(
  x = Idents(sce_cluster),
  from = names(new_cluster_ids),
  to = new_cluster_ids
)

celltype_colors <- c(
  "Schwann" = "#AEC7E8",
  "Fibroblast"= "#FFBB78",
  "Endothelial" = "#98DF8A",
  "Epithelial" = "#F7B6D2",
  "Stellate" = "#DBDB8D"
)

pdf(file.path(output_dir, "UMAP_annotation_PAAD_nonimmune.pdf"), width = 10, height = 8)
DimPlot(sce_cluster, group.by = "manual_annotation",label = TRUE,repel = TRUE,cols = celltype_colors) +
  ggtitle("PAAD Nonimmune Cell Annotation") 
dev.off()

saveRDS(sce_cluster, "my_data/PAAD_nonimmunecell_annotated_main.rds")



####################### 1.3 Convert barcode × gene matrix to celltype × gene matrix, aggregate expression values by cell type
expression_matrix <- GetAssayData(sce_cluster, assay = "RNA", slot = "counts")
expr_df <- as.data.frame(as.matrix(t(expression_matrix)))  
expr_df$celltype <- sce_cluster$manual_annotation

celltype_sum_matrix <- expr_df %>%
  group_by(celltype) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("celltype") %>%
  t()  

write.csv(celltype_sum_matrix, file= file.path(output_dir, "celltype_sum_countexpression_PAAD_nonimmune.csv"))





################################################## 2. Further Subpopulation Analysis
### 2.1 Epithelial Cell Subpopulation Analysis
Epithelial_cells <- subset(sce_cluster, subset = manual_annotation == "Epithelial")

# Re-preprocessing
Epithelial_cells <- NormalizeData(Epithelial_cells, ormalization.method = "LogNormalize", scale.factor = 10000)
Epithelial_cells <- FindVariableFeatures(Epithelial_cells, selection.method = "vst", nfeatures = 3000)
Epithelial_cells <- ScaleData(Epithelial_cells)
Epithelial_cells <- RunPCA(Epithelial_cells, features = VariableFeatures(Epithelial_cells))
ElbowPlot(Epithelial_cells, ndims = 50)
cumulative_var <- cumsum(Epithelial_cells@reductions$pca@stdev^2 / sum(Epithelial_cells@reductions$pca@stdev^2))
elbow_point <- which.max(diff(cumulative_var) < 0.01)  
print(paste("Estimated elbow point at PC:", elbow_point))
plot(cumulative_var, type = "b", xlab = "PCs", ylab = "Cumulative Variance Explained")
abline(v = 22, col = "blue", lty = 2)
abline(h = 0.8, col = "grey", lty = 3)  
Epithelial_cells <- FindNeighbors(Epithelial_cells, reduction = "pca", dims = 1:22)
Epithelial_cells <- RunUMAP(Epithelial_cells, reduction = "pca", dims = 1:22)

#Clustering
Epithelial_cells <- FindClusters(Epithelial_cells, resolution = 1.5, algorithm = 4, random.seed = 42)
table(Epithelial_cells@meta.data$seurat_clusters)

DimPlot(Epithelial_cells, reduction = "umap", label = TRUE) + NoLegend()

#Find marker genes for Epithelial cell subpopulations
Epithelial_markers <- FindAllMarkers(
  Epithelial_cells,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 1,
  test.use = "wilcox"
)
top10_Epithelialmarkers <- Epithelial_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(Epithelial_markers, file= file.path(output_dir, "Epithelial_cluster_markers_PAAD_nonimmune.csv"))
write.csv(top10_Epithelialmarkers, file= file.path(output_dir, "Top10Epithelial_cluster_markers_PAAD_nonimmune.csv"))


#Use the GO enrichment terms of each cluster to roughly determine the names of the annotations
for (cl in 1:23) {
  cluster_genes <- subset(Epithelial_markers, cluster == cl)$gene
  ego <- enrichGO(
    gene = cluster_genes,
    OrgDb = "org.Hs.eg.db",
    keyType = "SYMBOL",
    ont = "BP", 
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  if (nrow(ego) > 0) {
    pdf_file <- file.path(output_dir, paste0("Epithelial_GO_markers_cluster", cl, "_PAAD_nonimmune.pdf"))
    pdf(pdf_file, width = 12, height = 10)
    print(dotplot(ego))
    dev.off()
  } else {
    message("Cluster ", cl, " has no significant enrichment results, skipping.")
  }
}

#Create Epithelial cell marker sets
Epithelial_subset_markers <- list(
  EMTEp=c("VIM","ZEB1","ZEB2","SNAI2","FN1"),
  metEp=c("HK2","SLC2A1","VEGFA","SLC1A5","GLUD1","BCAT1"),
  imEp=c("HLA-DRA","HLA-DRB1","HLA-DPA1","CD74","CTSS","CXCL1","CCL2"),
  pEp=c("MKI67","TOP2A","CCNB1","CCNB2","PCNA","AURKA"),
  sEp=c("INSL4","CCK","MSMB","ECEL1","CGA","CST2","CST4","CST5","BPI","ALB","INSL4","MUC1","MUC3A","SPINK1","CTSE","LGALS4"),
  ciliumEp=c("C9orf24","CFAP126","C5orf49","DRC1"),
  sarcomatoidEp=c("COL27A1","MYBPC1")
  
)

FeaturePlot(Epithelial_cells, features = unlist(Epithelial_subset_markers), ncol = 4)
DotPlot(Epithelial_cells, features = unique(unlist(Epithelial_subset_markers))) + 
  RotatedAxis() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red")
VlnPlot(Epithelial_cells, features = unique(unlist(Epithelial_subset_markers)), stack = T, flip = T) + 
  NoLegend()

#Annotation
new_Epithelial_labels <- c(
  "1" = "metEp",
  "2" = "sEp",
  "3" = "imEp",
  "4" = "sEp",
  "5" = "EMTEp",
  "6" = "sEp",
  "7" = "sEp",
  "8" = "imEp",
  "9" = "imEp",
  "10" = "EMTEp",
  "11" = "pEp",
  "12" = "sarcomatoidEp",
  "13" = "imEp",
  "14" = "ciliumEp",
  "15" = "EMTEp",
  "16" = "metEp",
  "17" = "sEp",
  "18" = "sEp",
  "19" = "imEp",
  "20" = "imEp",
  "21" = "imEp",
  "22" = "metEp",
  "23" = "imEp"
)

Epithelial_cells$Epithelial_subset_annotation <- plyr::mapvalues(
  x = Idents(Epithelial_cells),
  from = names(new_Epithelial_labels),
  to = new_Epithelial_labels
)

celltype_colors <- c(
  "EMTEp" = "#AEC7E8",
  "metEp"= "#FFBB78",
  "imEp" = "#98DF8A",
  "pEp" = "#F7B6D2",
  "sEp" = "#DBDB8D",
  "ciliumEp"="#8C564B",
  "sarcomatoidEp"="#17BECF"
)

pdf(file.path(output_dir, "Epithelial_UMAP_annotation_PAAD_nonimmune.pdf"), width = 8, height = 6)
DimPlot(Epithelial_cells, group.by = "Epithelial_subset_annotation", label = TRUE, repel = TRUE,cols = celltype_colors) +
  ggtitle("Epithelial Subset Cell Annotation") 
dev.off()

saveRDS(Epithelial_cells, "my_data/PAAD_nonimmunecell_annotated_Epithelialsubset.rds")

#Convert Epithelial cell barcode × gene matrix to celltype × gene matrix for subsequent ferroptosis scoring
expression_matrix_Epithelial <- GetAssayData(Epithelial_cells, assay = "RNA", slot = "counts")
expr_df_Epithelial <- as.data.frame(as.matrix(t(expression_matrix_Epithelial)))  
expr_df_Epithelial$celltype <- Epithelial_cells$Epithelial_subset_annotation

celltype_sum_matrix_Epithelial <- expr_df_Epithelial %>%
  group_by(celltype) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("celltype") %>%
  t()  

write.csv(celltype_sum_matrix_Epithelial, file= file.path(output_dir, "celltype_sum_countexpression_Epithelialsubset_PAAD_nonimmune.csv"))



################2.2 Fibroblast Cell Subpopulation Analysis
Fibroblast_cells <- subset(sce_cluster, subset = manual_annotation == "Fibroblast")

Fibroblast_cells <- NormalizeData(Fibroblast_cells, ormalization.method = "LogNormalize", scale.factor = 10000)
Fibroblast_cells <- FindVariableFeatures(Fibroblast_cells, selection.method = "vst", nfeatures = 3000)
Fibroblast_cells <- ScaleData(Fibroblast_cells)
Fibroblast_cells <- RunPCA(Fibroblast_cells, features = VariableFeatures(Fibroblast_cells))
ElbowPlot(Fibroblast_cells, ndims = 50)
cumulative_var <- cumsum(Fibroblast_cells@reductions$pca@stdev^2 / sum(Fibroblast_cells@reductions$pca@stdev^2))
elbow_point <- which.max(diff(cumulative_var) < 0.01)  
print(paste("Estimated elbow point at PC:", elbow_point))
plot(cumulative_var, type = "b", xlab = "PCs", ylab = "Cumulative Variance Explained")
abline(v = 22, col = "blue", lty = 2)
abline(h = 0.8, col = "grey", lty = 3)  
Fibroblast_cells <- FindNeighbors(Fibroblast_cells, reduction = "pca", dims = 1:22)
Fibroblast_cells <- RunUMAP(Fibroblast_cells, reduction = "pca", dims = 1:22)


Fibroblast_cells <- FindClusters(Fibroblast_cells, resolution = 1.5, algorithm = 4, random.seed = 42)
table(Fibroblast_cells@meta.data$seurat_clusters)

DimPlot(Fibroblast_cells, reduction = "umap", label = TRUE) + NoLegend()


Fibroblast_markers <- FindAllMarkers(
  Fibroblast_cells,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 1,
  test.use = "wilcox"
)
top10_Fibroblastmarkers <- Fibroblast_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(Fibroblast_markers, file= file.path(output_dir, "Fibroblastcell_cluster_markers_PAAD_nonimmune.csv"))
write.csv(top10_Fibroblastmarkers, file= file.path(output_dir, "Top10Fibroblastcell_cluster_markers_PAAD_nonimmune.csv"))


for (cl in 1:23) {
  cluster_genes <- subset(Fibroblast_markers, cluster == cl)$gene
  ego <- enrichGO(
    gene = cluster_genes,
    OrgDb = "org.Hs.eg.db",
    keyType = "SYMBOL",
    ont = "BP", 
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  if (nrow(ego) > 0) {
    pdf_file <- file.path(output_dir, paste0("Fibroblast_GO_markers_cluster", cl, "_PAAD_nonimmune.pdf"))
    pdf(pdf_file, width = 12, height = 10)
    print(dotplot(ego))
    dev.off()
  } else {
    message("Cluster ", cl, " has no significant enrichment results, skipping.")
  }
}


Fibroblast_subset_markers <- list(
  mCAF=c("FN1","LUM","DCN","VCAN","FAP","MMP11","LOXL2"),
  iCAF=c("IL6","PDGFRA","LIF","CXCL12","CCL2","CXCL1","CXCL2","HAS1"),
  hspCAF=c("HSPA6","HSP90AA1","HSPB1","DNAJB1","HMOX1","CRYAB","XBP1","ATF4"),
  dCAF=c("FRZB","SOX9","WT1","RUNX2"),
  nCAF=c("SBSPON","CAPN6","KLK1","SLCO1C1")
)

FeaturePlot(Fibroblast_cells, features = unlist(Fibroblast_subset_markers), ncol = 4)
DotPlot(Fibroblast_cells, features = unique(unlist(Fibroblast_subset_markers))) + 
  RotatedAxis() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red")
VlnPlot(Fibroblast_cells, features = unique(unlist(Fibroblast_subset_markers)), stack = T, flip = T) + 
  NoLegend()


new_Fibroblast_labels <- c(
  "1" = "mCAF",
  "2" = "dCAF",
  "3" = "mCAF",
  "4" = "mCAF",
  "5" = "dCAF",
  "6" = "mCAF",
  "7" = "dCAF",
  "8" = "mCAF",
  "9" = "iCAF",
  "10" = "mCAF",
  "11" = "mCAF",
  "12" = "dCAF",
  "13" = "dCAF",
  "14" = "hspCAF",
  "15" = "mCAF",
  "16" = "dCAF",
  "17" = "dCAF",
  "18" = "mCAF",
  "19" = "nCAF",
  "20" = "iCAF",
  "21" = "iCAF",
  "22" = "mCAF",
  "23" = "dCAF"
  
)

Fibroblast_cells$Fibroblast_subset_annotation <- plyr::mapvalues(
  x = Idents(Fibroblast_cells),
  from = names(new_Fibroblast_labels),
  to = new_Fibroblast_labels
)

celltype_colors <- c(
  "mCAF" = "#AEC7E8",
  "iCAF"= "#FFBB78",
  "hspCAF" = "#98DF8A",
  "dCAF" = "#F7B6D2",
  "nCAF" = "#DBDB8D"
)

pdf(file.path(output_dir, "Fibroblast_UMAP_annotation_PAAD_nonimmune.pdf"), width = 8, height = 6)
DimPlot(Fibroblast_cells, group.by = "Fibroblast_subset_annotation", label = TRUE, repel = TRUE,cols = celltype_colors) +
  ggtitle("Fibroblast Subset Cell Annotation")  
dev.off()

saveRDS(Fibroblast_cells, "my_data/PAAD_nonimmunecell_annotated_Fibroblastsubset.rds")


expression_matrix_Fibroblast <- GetAssayData(Fibroblast_cells, assay = "RNA", slot = "counts")
expr_df_Fibroblast <- as.data.frame(as.matrix(t(expression_matrix_Fibroblast)))  
expr_df_Fibroblast$celltype <- Fibroblast_cells$Fibroblast_subset_annotation

celltype_sum_matrix_Fibroblast <- expr_df_Fibroblast %>%
  group_by(celltype) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("celltype") %>%
  t()  

write.csv(celltype_sum_matrix_Fibroblast, file= file.path(output_dir, "celltype_sum_countexpression_Fibroblastsubset_PAAD_nonimmune.csv"))



################2.3 Stellate Cell Subpopulation Analysis
Stellate_cells <- subset(sce_cluster, subset = manual_annotation == "Stellate")


Stellate_cells <- NormalizeData(Stellate_cells, ormalization.method = "LogNormalize", scale.factor = 10000)
Stellate_cells <- FindVariableFeatures(Stellate_cells, selection.method = "vst", nfeatures = 3000)
Stellate_cells <- ScaleData(Stellate_cells)
Stellate_cells <- RunPCA(Stellate_cells, features = VariableFeatures(Stellate_cells))
ElbowPlot(Stellate_cells, ndims = 50)
cumulative_var <- cumsum(Stellate_cells@reductions$pca@stdev^2 / sum(Stellate_cells@reductions$pca@stdev^2))
elbow_point <- which.max(diff(cumulative_var) < 0.01)  
print(paste("Estimated elbow point at PC:", elbow_point))
plot(cumulative_var, type = "b", xlab = "PCs", ylab = "Cumulative Variance Explained")
abline(v = 27, col = "blue", lty = 2)
abline(h = 0.8, col = "grey", lty = 3)  
Stellate_cells <- FindNeighbors(Stellate_cells, reduction = "pca", dims = 1:27)
Stellate_cells <- RunUMAP(Stellate_cells, reduction = "pca", dims = 1:27)


Stellate_cells <- FindClusters(Stellate_cells, resolution = 2, algorithm = 4, random.seed = 42)
table(Stellate_cells@meta.data$seurat_clusters)

DimPlot(Stellate_cells, reduction = "umap", label = TRUE) + NoLegend()


Stellate_markers <- FindAllMarkers(
  Stellate_cells,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 1,
  test.use = "wilcox"
)
top10_Stellatemarkers <- Stellate_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(Stellate_markers, file= file.path(output_dir, "Stellate_cluster_markers_PAAD_nonimmune.csv"))
write.csv(top10_Stellatemarkers, file= file.path(output_dir, "Top10Stellate_cluster_markers_PAAD_nonimmune.csv"))


for (cl in 1:16) {
  cluster_genes <- subset(Stellate_markers, cluster == cl)$gene
  ego <- enrichGO(
    gene = cluster_genes,
    OrgDb = "org.Hs.eg.db",
    keyType = "SYMBOL",
    ont = "BP", 
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  if (nrow(ego) > 0) {
    pdf_file <- file.path(output_dir, paste0("Stellate_GO_markers_cluster", cl, "_PAAD_nonimmune.pdf"))
    pdf(pdf_file, width = 12, height = 10)
    print(dotplot(ego))
    dev.off()
  } else {
    message("Cluster ", cl, " has no significant enrichment results, skipping.")
  }
}


Stellate_subset_markers <- list(
  fPSC=c("ACTA2","COL1A1","COL3A1","COL5A1","TGFB1","TIMP1","TIMP2","FN1","SMAD2","SMAD3","PDGFRA","PDGFRB","MMP2","MMP9"),
  metPSC=c("HK2","LDHA","GLS","GLUD1","SLC1A5","FASN","ACACA","CPT1A","PLIN2","HIF1A","MYC","PPARG","SREBF1","PPARGC1A"),
  iPSC=c("IL6","CCL2","CXCL12","CLU","IER3"),
  imPSC=c("CD274","CD80","CD86","CTSL","PSMB8","PSMB9","CD74","CIITA"),
  nPSC=c("NGFR","GFAP","EDN1","S100B","NGF","NOS1","GDNF","PTGER3","PTGER1","SCN4B"),
  cPSC=c("MYH11","TAGLN","CNN1","ACTG2","FLNA","LOX","ITGB1")
)

FeaturePlot(Stellate_cells, features = unlist(Stellate_subset_markers), ncol = 4)
DotPlot(Stellate_cells, features = unique(unlist(Stellate_subset_markers))) + 
  RotatedAxis() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red")
VlnPlot(Stellate_cells, features = unique(unlist(Stellate_subset_markers)), stack = T, flip = T) + 
  NoLegend()


new_Stellate_labels <- c(
  "1" = "nPSC",
  "2" = "fPSC",
  "3" = "cPSC",
  "4" = "cPSC",
  "5" = "iPSC",
  "6" = "metPSC",
  "7" = "imPSC",
  "8" = "iPSC",
  "9" = "imPSC",
  "10" = "iPSC",
  "11" = "imPSC",
  "12" = "cPSC",
  "13" = "iPSC",
  "14" = "imPSC",
  "15" = "fPSC",
  "16" = "imPSC"
)

Stellate_cells$Stellate_subset_annotation <- plyr::mapvalues(
  x = Idents(Stellate_cells),
  from = names(new_Stellate_labels),
  to = new_Stellate_labels
)

celltype_colors <- c(
  "fPSC" = "#AEC7E8",
  "metPSC"= "#FFBB78",
  "iPSC" = "#98DF8A",
  "imPSC" = "#F7B6D2",
  "nPSC" = "#DBDB8D",
  "cPSC" = "#8C564B"
)

pdf(file.path(output_dir, "Stellate_UMAP_annotation_PAAD_nonimmune.pdf"), width = 8, height = 6)
DimPlot(Stellate_cells, group.by = "Stellate_subset_annotation", label = TRUE, repel = TRUE,cols = celltype_colors) +
  ggtitle("Stellate Subset Cell Annotation")  
dev.off()

saveRDS(Stellate_cells, "my_data/PAAD_nonimmunecell_annotated_Stellatesubset.rds")


expression_matrix_Stellate <- GetAssayData(Stellate_cells, assay = "RNA", slot = "counts")
expr_df_Stellate <- as.data.frame(as.matrix(t(expression_matrix_Stellate)))  
expr_df_Stellate$celltype <- Stellate_cells$Stellate_subset_annotation

celltype_sum_matrix_Stellate <- expr_df_Stellate %>%
  group_by(celltype) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("celltype") %>%
  t()  

write.csv(celltype_sum_matrix_Stellate, file= file.path(output_dir, "celltype_sum_countexpression_Stellatesubset_PAAD_nonimmune.csv"))



################2.4 Endothelial Cell Subpopulation Analysis
Endothelial_cells <- subset(sce_cluster, subset = manual_annotation == "Endothelial")

Endothelial_cells <- NormalizeData(Endothelial_cells, ormalization.method = "LogNormalize", scale.factor = 10000)
Endothelial_cells <- FindVariableFeatures(Endothelial_cells, selection.method = "vst", nfeatures = 3000)
Endothelial_cells <- ScaleData(Endothelial_cells)
Endothelial_cells <- RunPCA(Endothelial_cells, features = VariableFeatures(Endothelial_cells))
ElbowPlot(Endothelial_cells, ndims = 50)
cumulative_var <- cumsum(Endothelial_cells@reductions$pca@stdev^2 / sum(Endothelial_cells@reductions$pca@stdev^2))
elbow_point <- which.max(diff(cumulative_var) < 0.01)  
print(paste("Estimated elbow point at PC:", elbow_point))
plot(cumulative_var, type = "b", xlab = "PCs", ylab = "Cumulative Variance Explained")
abline(v = 25, col = "blue", lty = 2)
abline(h = 0.8, col = "grey", lty = 3)  
Endothelial_cells <- FindNeighbors(Endothelial_cells, reduction = "pca", dims = 1:25)
Endothelial_cells <- RunUMAP(Endothelial_cells, reduction = "pca", dims = 1:25)


Endothelial_cells <- FindClusters(Endothelial_cells, resolution = 2, algorithm = 4, random.seed = 42)
table(Endothelial_cells@meta.data$seurat_clusters)

DimPlot(Endothelial_cells, reduction = "umap", label = TRUE) + NoLegend()


Endothelial_markers <- FindAllMarkers(
  Endothelial_cells,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 1,
  test.use = "wilcox"
)
top10_Endothelialmarkers <- Endothelial_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.csv(Endothelial_markers, file= file.path(output_dir, "Endothelial_cluster_markers_PAAD_nonimmune.csv"))
write.csv(top10_Endothelialmarkers, file= file.path(output_dir, "Top10Endothelial_cluster_markers_PAAD_nonimmune.csv"))


for (cl in 1:19) {
  cluster_genes <- subset(Endothelial_markers, cluster == cl)$gene
  ego <- enrichGO(
    gene = cluster_genes,
    OrgDb = "org.Hs.eg.db",
    keyType = "SYMBOL",
    ont = "BP", 
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  if (nrow(ego) > 0) {
    pdf_file <- file.path(output_dir, paste0("Endothelial_GO_markers_cluster", cl, "_PAAD_nonimmune.pdf"))
    pdf(pdf_file, width = 12, height = 10)
    print(dotplot(ego))
    dev.off()
  } else {
    message("Cluster ", cl, " has no significant enrichment results, skipping.")
  }
}


Endothelial_subset_markers <- list(
  lEn=c("LYVE1","PROX1","PDPN","FLT4"),
  vEn_Artery=c("CD34","VWF","CLDN5","EFNB2","GJA5","DLL4","SOX17"),
  vEn_Vein=c("CD34","VWF","CLDN5","NR2F2","ACKR1","EPHB4"),
  vEn_Capillary=c("CD34","VWF","CLDN5","RGS5","PLVAP","ESM1"),
  iEn=c("ANGPT2","SELE","VCAM1","VEGFA","SLC2A1","XBP1","ICAM1"),
  EMTEn=c("CD34","ACTA2","FAP","SNAI1","SNAI2","COL1A1")
  
)

FeaturePlot(Endothelial_cells, features = unlist(Endothelial_subset_markers), ncol = 4)
DotPlot(Endothelial_cells, features = unique(unlist(Endothelial_subset_markers))) + 
  RotatedAxis() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red")
VlnPlot(Endothelial_cells, features = unique(unlist(Endothelial_subset_markers)), stack = T, flip = T) + 
  NoLegend()


new_Endothelial_labels <- c(
  "1" = "vEn_Vein",
  "2" = "iEn",
  "3" = "iEn",
  "4" = "iEn",
  "5" = "vEn_Vein",
  "6" = "vEn_Artery",
  "7" = "vEn_Artery",
  "8" = "vEn_Artery",
  "9" = "vEn_Artery",
  "10" = "vEn_Artery",
  "11" = "vEn_Capillary",
  "12" = "vEn_Artery",
  "13" = "EMTEn",
  "14" = "iEn",
  "15" = "vEn_Artery",
  "16" = "vEn_Vein",
  "17" = "EMTEn",
  "18" = "lEn",
  "19" = "vEn_Capillary"
  
)

Endothelial_cells$Endothelial_subset_annotation <- plyr::mapvalues(
  x = Idents(Endothelial_cells),
  from = names(new_Endothelial_labels),
  to = new_Endothelial_labels
)

celltype_colors <- c(
  "lEn" = "#AEC7E8",
  "vEn_Artery"= "#FFBB78",
  "vEn_Vein" = "#98DF8A",
  "vEn_Capillary" = "#F7B6D2",
  "iEn" = "#DBDB8D",
  "EMTEn" = "#8C564B"
)

pdf(file.path(output_dir, "Endothelial_UMAP_annotation_PAAD_nonimmune.pdf"), width = 8, height = 6)
DimPlot(Endothelial_cells, group.by = "Endothelial_subset_annotation", label = TRUE, repel = TRUE,cols = celltype_colors) +
  ggtitle("Endothelial Subset Cell Annotation")   
dev.off()

saveRDS(Endothelial_cells, "my_data/PAAD_nonimmunecell_annotated_Endothelialsubset.rds")


expression_matrix_Endothelial <- GetAssayData(Endothelial_cells, assay = "RNA", slot = "counts")
expr_df_Endothelial <- as.data.frame(as.matrix(t(expression_matrix_Endothelial)))  
expr_df_Endothelial$celltype <- Endothelial_cells$Endothelial_subset_annotation

celltype_sum_matrix_Endothelial <- expr_df_Endothelial %>%
  group_by(celltype) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("celltype") %>%
  t()  

write.csv(celltype_sum_matrix_Endothelial, file= file.path(output_dir, "celltype_sum_countexpression_Endothelialsubset_PAAD_nonimmune.csv"))





