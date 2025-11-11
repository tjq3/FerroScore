########scRNA-seq data preprocess
#############################################



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Seurat","Matrix","ggplot2","scDblFinder"))
library(Seurat)
library(Matrix)
library(ggplot2)
library(scDblFinder)




#Load barcodes, features, .mtx files
k <- Matrix::readMM("my_data/matrix.mtx")
mtx <- as(k, "CsparseMatrix")
dim(mtx)

cl <- read.table('my_data/barcodes.tsv',header = FALSE, stringsAsFactors = FALSE)[,1]
rl <- read.table('my_data/features.tsv',header = FALSE, stringsAsFactors = FALSE)[,2]

#Convert transcript IDs like AL627309.1, AL627309.2, AL627309.3 in features to gene level AL627309, and take sum as expression
gene_names <- sub("\\.[0-9]+$", "", rl)
unique_genes <- unique(gene_names)
agg_mtx <- matrix(0, nrow = length(unique_genes), ncol = ncol(mtx))

for(i in seq_along(unique_genes)) {
  idx <- which(gene_names == unique_genes[i])
  agg_mtx[i,] <- colSums(mtx[idx,, drop = FALSE])
}

#Organize dataset
rownames(agg_mtx) <- unique_genes
colnames(agg_mtx) <- cl
sce <- CreateSeuratObject(counts = as(agg_mtx, "dgCMatrix"))

saveRDS(sce, "my_data/PAAD_immunecell.rds")




#########Data Preprocessing
#step1: Quality Control Filtering
sce[["percent.mt"]] <- PercentageFeatureSet(sce, assay = "RNA", pattern = "^MT-")

sce_sce <- as.SingleCellExperiment(sce)
set.seed(42)
sce_sce <- scDblFinder(sce_sce)
sce[["scDblFinder_score"]] <- sce_sce$scDblFinder.score
sce[["scDblFinder_class"]] <- sce_sce$scDblFinder.class

sce_qc_single <- subset(sce, 
                        subset = nCount_RNA > 2000 &       #UMI > 2000
                          nFeature_RNA > 500 &             #Number of genes > 500
                          nFeature_RNA < 7000 &            #Number of genes < 7000
                          percent.mt < 10 &                #Mitochondrial genes < 10%
                          scDblFinder_class == "singlet")  #Keep only single cells

print(paste("Original number of cells:", ncol(sce)))
print(paste("Number of cells after filtering:", ncol(sce_qc_single)))
dim(sce_qc_single)
head(sce_qc_single@meta.data)


#step2: Data Normalization
sce_lognormal_single <- NormalizeData(sce_qc_single, 
                                      normalization.method = "LogNormalize",
                                      scale.factor = 10000)
dim(sce_lognormal_single)


#step3: Feature Selection - Highly Variable Genes
sce_vst_single <- FindVariableFeatures(sce_lognormal_single, 
                                       selection.method = "vst",
                                       nfeatures = 3000)
dim(sce_vst_single)
top20_single <- head(VariableFeatures(sce_vst_single), 20)
length(VariableFeatures(sce_vst_single))


#step4: Data Scaling
sce_scale_single <- ScaleData(sce_vst_single)
dim(sce_scale_single)


#step5: Linear Dimensionality Reduction - PCA
sce_pca_single <- RunPCA(sce_scale_single, 
                         features = VariableFeatures(sce_scale_single))
dim(sce_pca_single)
ElbowPlot(sce_pca_single, ndims = 50)

cumulative_var_single <- cumsum(sce_pca_single@reductions$pca@stdev^2 / sum(sce_pca_single@reductions$pca@stdev^2))
elbow_point_single <- which.max(diff(cumulative_var_single) < 0.01)  
print(paste("Estimated elbow point at PC:", elbow_point_single))
plot(cumulative_var_single, type = "b", xlab = "PCs", ylab = "Cumulative Variance Explained")
abline(v = elbow_point, col = "blue", lty = 2)
abline(h = 0.8, col = "grey", lty = 3)  
ncol(sce_pca_single@reductions$pca@cell.embeddings)


#step6: Cell Clustering
sce_nei_single <- FindNeighbors(sce_pca_single, dims = 1:21)  
dim(sce_nei_single)


#step7: Nonlinear Dimensionality Reduction - UMAP
sce_umap_single <- RunUMAP(sce_nei_single, dims = 1:21)
dim(sce_umap_single)
DimPlot(sce_umap_single, reduction = "umap")


saveRDS(sce_umap_single, "my_data/PAAD_immunecell_datapre.rds")








