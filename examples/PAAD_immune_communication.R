######Ferroptosis drives an immunosuppressive microenvironment through intercellular communication
###################################################################



usethis::create_github_token()
usethis::edit_r_environ()
Sys.getenv("GITHUB_TOKEN")

install.packages("devtools")
devtools::install_github("jinworks/CellChat")
devtools::install_github('immunogenomics/presto')
install.packages("reticulate")

library(Seurat)
library(CellChat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(patchwork)
options(stringsAsFactors = FALSE)
library(reticulate)
py_config()
py_install("umap-learn")
library(igraph)
library(readxl)



output_dir <- "my_results/immunecommunication"




#######Data Loading
immune <- readRDS("my_data/PAAD_immunecell_annotated_detailed.rds")

nonimmune_all <- readRDS("my_data/PAAD_nonimmunecell_annotated_main.rds")
nonimmune_sch <- subset(nonimmune_all, subset = manual_annotation == "Schwann")
saveRDS(nonimmune_sch, "my_data/PAAD_nonimmunecell_Schwannsubset.rds")
nonimmune_Ep<-readRDS("my_data/PAAD_nonimmunecell_annotated_Epithelialsubset.rds")
nonimmune_CAF<-readRDS("my_data/PAAD_nonimmunecell_annotated_Fibroblastsubset.rds")
nonimmune_PSC<-readRDS("my_data/PAAD_nonimmunecell_annotated_Stellatesubset.rds")
nonimmune_En<-readRDS("my_data/PAAD_nonimmunecell_annotated_Endothelialsubset.rds")


immune@meta.data <- immune@meta.data %>%
  mutate(
    cell_origin = "immune",
    major_celltype = manual_annotation,
    subset_annotation = detailed_annotation
  )
nonimmune_sch@meta.data <- nonimmune_sch@meta.data %>%
  mutate(
    cell_origin = "nonimmune",
    major_celltype = manual_annotation,
    subset_annotation = manual_annotation
  )
nonimmune_Ep@meta.data <- nonimmune_Ep@meta.data %>%
  mutate(
    cell_origin = "nonimmune",
    major_celltype = manual_annotation,
    subset_annotation = Epithelial_subset_annotation
  )
nonimmune_CAF@meta.data <- nonimmune_CAF@meta.data %>%
  mutate(
    cell_origin = "nonimmune",
    major_celltype = manual_annotation,
    subset_annotation = Fibroblast_subset_annotation
  )
nonimmune_PSC@meta.data <- nonimmune_PSC@meta.data %>%
  mutate(
    cell_origin = "nonimmune",
    major_celltype = manual_annotation,
    subset_annotation = Stellate_subset_annotation
  )
nonimmune_En@meta.data <- nonimmune_En@meta.data %>%
  mutate(
    cell_origin = "nonimmune",
    major_celltype = manual_annotation,
    subset_annotation = Endothelial_subset_annotation
  )

##Unify meta information
add_full_annotation <- function(obj) {
  meta <- obj@meta.data
  if ("detailed_annotation" %in% colnames(meta)) {
    meta$full_annotation <- meta$detailed_annotation #immune cell
  } else if ("Epithelial_subset_annotation" %in% colnames(meta)) {
    meta$full_annotation <- paste("Ep", meta$Epithelial_subset_annotation, sep = "_") #Ep cell
  } else if ("Fibroblast_subset_annotation" %in% colnames(meta)) {
    meta$full_annotation <- paste("CAF", meta$Fibroblast_subset_annotation, sep = "_") #CAF cell
  } else if ("Stellate_subset_annotation" %in% colnames(meta)) {
    meta$full_annotation <- paste("PSC", meta$Stellate_subset_annotation, sep = "_") #PSC cell
  } else if ("Endothelial_subset_annotation" %in% colnames(meta)) {
    meta$full_annotation <- paste("En", meta$Endothelial_subset_annotation, sep = "_") #En cell
  } else {
    meta$full_annotation <- meta$manual_annotation #sch cell
  }
  obj@meta.data <- meta
  return(obj)
}

immune <- add_full_annotation(immune)
nonimmune_sch <- add_full_annotation(nonimmune_sch)
nonimmune_Ep <- add_full_annotation(nonimmune_Ep)
nonimmune_CAF <- add_full_annotation(nonimmune_CAF)
nonimmune_PSC <- add_full_annotation(nonimmune_PSC)
nonimmune_En <- add_full_annotation(nonimmune_En)

#Unify gene
common_genes <- Reduce(intersect, list(
  rownames(immune),
  rownames(nonimmune_sch), 
  rownames(nonimmune_Ep),
  rownames(nonimmune_CAF),
  rownames(nonimmune_PSC),
  rownames(nonimmune_En)
))

immune <- immune[common_genes, ]
nonimmune_sch <- nonimmune_sch[common_genes, ] 
nonimmune_Ep <- nonimmune_Ep[common_genes, ]
nonimmune_CAF <- nonimmune_CAF[common_genes, ]
nonimmune_PSC <- nonimmune_PSC[common_genes, ]
nonimmune_En <- nonimmune_En[common_genes, ]

##Merge all datasets
combined <- merge(immune, y = list(nonimmune_sch, nonimmune_Ep, nonimmune_CAF, nonimmune_PSC, nonimmune_En))
combined@meta.data <- combined@meta.data %>%
  mutate(
    dataset_origin = case_when(
      grepl("^data1_", colnames(combined)) ~ "nonimmune",
      TRUE ~ "immune"
    )
  )
table(combined$full_annotation)

saveRDS(combined, "my_data/PAAD_Immunitycombined_dataset.rds")
 



#########Build global cell communication network
combined <- JoinLayers(combined, assay = "RNA")
data.input <- GetAssayData(combined, assay = "RNA", layer = "data")
cellchat <- createCellChat(object = data.input, meta = combined@meta.data, group.by = "full_annotation")

#Set ligand-receptor database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 
cellchat@DB <- CellChatDB.use

#Preprocessing to filter low-expressed genes
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#Calculate communication probability
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

#Integrate pathway-level analysis
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

#Visualization
pdf(file.path(output_dir, "global_circle.pdf"), width = 14, height = 8)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
dev.off()


pdf(file.path(output_dir, "single_circle.pdf"), width = 14, height = 40)
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()




#######Cell communication among high ferroptosis cells
##Communication network among high ferroptosis cells
high_ferro_cells <- c("Macrophage", "CD8 Tcm", "NK", "CAF_nCAF", "PSC_iPSC", "En_EMTEn", "Ep_imEp")
combined_ferro <- subset(combined, subset = full_annotation %in% high_ferro_cells)
table(combined_ferro$full_annotation)

combined_ferro <- JoinLayers(combined_ferro, assay = "RNA")
data.input_ferro <- GetAssayData(combined_ferro, assay = "RNA", layer = "data")
cellchat_ferro <- createCellChat(object = data.input_ferro, 
                                 meta = combined_ferro@meta.data, 
                                 group.by = "full_annotation")

CellChatDB <- CellChatDB.human 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 
cellchat_ferro@DB <- CellChatDB.use
cellchat_ferro <- subsetData(cellchat_ferro)
cellchat_ferro <- identifyOverExpressedGenes(cellchat_ferro)
cellchat_ferro <- identifyOverExpressedInteractions(cellchat_ferro)
cellchat_ferro <- computeCommunProb(cellchat_ferro, type = "triMean")
cellchat_ferro <- filterCommunication(cellchat_ferro, min.cells = 10)
cellchat_ferro <- computeCommunProbPathway(cellchat_ferro)
cellchat_ferro <- aggregateNet(cellchat_ferro)

pdf(file.path(output_dir, "highferrocell_global_circle.pdf"), width = 14, height = 8)
groupSize <- as.numeric(table(cellchat_ferro@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_ferro@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_ferro@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
dev.off()
pdf(file.path(output_dir, "highferrocell_single_circle.pdf"), width = 14, height = 40)
mat <- cellchat_ferro@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()


##Generated signaling pathways
cellchat_ferro@netP$pathways
write.csv(as.data.frame(cellchat_ferro@netP$pathways), file= file.path(output_dir, "signaling_allsignaling.csv"))

##Identify which signals contribute most to incoming/outgoing communication in high ferroptosis cell populations
pdf(file.path(output_dir, "signaling_roles.pdf"), width = 14, height = 9)
cellchat_ferro <- netAnalysis_computeCentrality(cellchat_ferro, slot.name = "netP") 
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_ferro, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_ferro, pattern = "incoming")
ht1 + ht2
dev.off()

##Further extract key research pathways, identify sending/receiving/mediating/influencing cells in each signaling pathway
pdf(file.path(output_dir, "signaling_roles_heatmap.pdf"), width = 14, height = 8)
pathways.show <- c("MIF","SPP1","GALECTIN","MK","PLAU") 
netAnalysis_signalingRole_network(cellchat_ferro, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
dev.off()

##Examine cell communication among high ferroptosis cells for individual signaling pathways
pathways <- c("MIF", "SPP1", "GALECTIN", "MK", "PLAU")
for (p in pathways) {
  pdf(file.path(output_dir, paste0("signaling_", p, "-signaling_circle.pdf")), width = 10, height = 8)
  netVisual_aggregate(cellchat_ferro, signaling = p, layout = "circle")
  dev.off()
}

##Examine gene expression of individual signaling pathways in high ferroptosis cells
for (p in pathways) {
  pdf(file.path(output_dir, paste0("signaling_", p, "-signaling_geneexpression.pdf")), width = 10, height = 8)
  plotGeneExpression(cellchat_ferro, signaling = p, enriched.only = TRUE, type = "violin")
  dev.off()
}

##Contribution of each ligand-receptor pair in individual signaling pathways
for (p in pathways) {
  pdf(file.path(output_dir, paste0("LR_", p, "-signaling_LRcontribution.pdf")), width = 10, height = 6)
  netAnalysis_contribution(cellchat_ferro, signaling = p)
  dev.off()
}




###Ferroptosis "Cell-Gene-Signaling" autocrine/paracrine regulatory network
#Get ferroptosis score data
ferro_scores <- data.frame(
  cell_type = c("Macrophage", "CD8 Tcm", "NK", "CAF_nCAF", "PSC_iPSC", "En_EMTEn", "Ep_imEp"),
  score = c(1.436645, 1.081676, 0.7276803, 2.073027, 0.8332458, 0.6732866, 0.6421859)
)

#Get DEGs data related to ferroptosis
ferr_db_genes <- read.csv("my_data/ferroptosis_all.csv")
ferro_genes <- unique(ferr_db_genes$symbol)

diff_genes_Macro <- read_excel("my_data/diff_genes_PAAD_immune_detailed_Macro.xlsx")
sig_genes_Macro <- diff_genes_Macro[diff_genes_Macro$PValue < 0.05 & abs(diff_genes_Macro$logFC) > 1, ]
ferro_diff_genes_Macro <- sig_genes_Macro[sig_genes_Macro$Gene %in% ferro_genes, ]

diff_genes_CD8Tcm <- read_excel("my_data/diff_genes_PAAD_immune_detailed_CD8Tcm.xlsx")
sig_genes_CD8Tcm <- diff_genes_CD8Tcm[diff_genes_CD8Tcm$PValue < 0.05 & abs(diff_genes_CD8Tcm$logFC) > 1, ]
ferro_diff_genes_CD8Tcm <- sig_genes_CD8Tcm[sig_genes_CD8Tcm$Gene %in% ferro_genes, ]

diff_genes_NK <- read_excel("my_data/diff_genes_PAAD_immune_detailed_NK.xlsx")
sig_genes_NK <- diff_genes_NK[diff_genes_NK$PValue < 0.05 & abs(diff_genes_NK$logFC) > 1, ]
ferro_diff_genes_NK <- sig_genes_NK[sig_genes_NK$Gene %in% ferro_genes, ]

diff_genes_EMTEn <- read_excel("my_data/diff_genes_PAAD_nonimmune_Endothelial_EMTEn.xlsx")
sig_genes_EMTEn <- diff_genes_EMTEn[diff_genes_EMTEn$PValue < 0.05 & abs(diff_genes_EMTEn$logFC) > 1, ]
ferro_diff_genes_EMTEn <- sig_genes_EMTEn[sig_genes_EMTEn$Gene %in% ferro_genes, ]

diff_genes_imEp <- read_excel("my_data/diff_genes_PAAD_nonimmune_Epithelial_imEp.xlsx")
sig_genes_imEp <- diff_genes_imEp[diff_genes_imEp$PValue < 0.05 & abs(diff_genes_imEp$logFC) > 1, ]
ferro_diff_genes_imEp <- sig_genes_imEp[sig_genes_imEp$Gene %in% ferro_genes, ]

diff_genes_nCAF <- read_excel("my_data/diff_genes_PAAD_nonimmune_Fibroblast_nCAF.xlsx")
sig_genes_nCAF <- diff_genes_nCAF[diff_genes_nCAF$PValue < 0.05 & abs(diff_genes_nCAF$logFC) > 1, ]
ferro_diff_genes_nCAF <- sig_genes_nCAF[sig_genes_nCAF$Gene %in% ferro_genes, ]

diff_genes_iPSC <- read_excel("my_data/diff_genes_PAAD_nonimmune_Stellate_iPSC.xlsx")
sig_genes_iPSC <- diff_genes_iPSC[diff_genes_iPSC$PValue < 0.05 & abs(diff_genes_iPSC$logFC) > 1, ]
ferro_diff_genes_iPSC <- sig_genes_iPSC[sig_genes_iPSC$Gene %in% ferro_genes, ]

ferro_diff_genes <- rbind(
  ferro_diff_genes_Macro,
  ferro_diff_genes_CD8Tcm,
  ferro_diff_genes_NK,
  ferro_diff_genes_EMTEn,
  ferro_diff_genes_imEp,
  ferro_diff_genes_nCAF,
  ferro_diff_genes_iPSC
)
ferro_diff_genes <- unique(ferro_diff_genes)




###### Autocrine regulatory networks for each cell type
target_cell_types <- c("Macrophage", "CD8_Tcm", "NK", "CAF_nCAF", 
                       "PSC_iPSC", "En_EMTEn", "Ep_imEp")
actual_col_names <- c("Macrophage", "CD8 Tcm", "NK", "CAF-nCAF",
                      "PSC-iPSC", "En-EMTEn", "Ep-imEp")

diff_genes_list <- list(
  Macrophage = ferro_diff_genes_Macro$Gene,
  CD8_Tcm = ferro_diff_genes_CD8Tcm$Gene,
  NK = ferro_diff_genes_NK$Gene,
  CAF_nCAF = ferro_diff_genes_nCAF$Gene,
  PSC_iPSC = ferro_diff_genes_iPSC$Gene,
  En_EMTEn = ferro_diff_genes_EMTEn$Gene,
  Ep_imEp = ferro_diff_genes_imEp$Gene
)

gene_pathway_mappings <- list(
  Macrophage = list(
    "MT1G" = c("TNF"),#all
    "ABHD12" = c("EDN","EGF","VEGF","GAS","IGF","PDGF","CXCL","PROS"), #all
    "FABP4" = c("CXCL","IGF","IGFBP","NGF","BMP","FGF","PDGF","SPP1","EGF","IL6","BAG","CCL","PLAU","IFNG","ANGPT","VEGF","CSF","EDN"),#all
    "FTL" = c("GRN","CXCL","EGF","PDGF","KLK","TNF"),#all
    "CTSB" = c("BMP", "KLK","PROS","GRN","CXCL","IL6","FASLG","IFNG","CCL","TNF","IGF","NGF","EGF"),#all
    "LGMN" = c("GRN", "CXCL","EGF","IGF","IGFBP"),#all
    "LAMP2" = c("EGF", "VEGF","IGF","PDGF","PROS","GRN","IL6","TNF","IFNG","BMP","SEMA3A","CXCL","PLAU","IGFBP","CSF"),#all
    "RB1" = c("EGF","FGF","VEGF","IGF","FASLG","CSF"),#top20
    "SLC38A1" = c("IGF", "IGFBP","TNF","IL6","IFNG"),#all
    "HMOX1" = c("EGF", "PDGF","CCL","FGF","CXCL","IFNG","IL16","GRN","IL6","TNF","VEGF","CSF","MIF","KLK","EDN"),#top20
    "PRDX1" = c("EGF", "PDGF","FASLG","KLK","TNF","CCL","IGF","IGFBP","CSF"),#all
    "PLIN2" = c("ANGPTL","TNF","IGFBP","CXCL","BMP","IL6","IGF","CSF","EDN","MIF","OSM")#all
  ),
  CD8_Tcm = list(
    "BRD2" = c("BMP","GDF","CSF","BAG"),#all
    "SLC1A5" = c("IGF","IGFBP","CXCL","CCL","SPP1","CSF")#all
  ),
  NK = list(
    "TGFB1" = c("BMP","CXCL","FGF","VEGF","GAS","GDF","IGF","IL6","NGF","PDGF","SEMA3","EGF","FASLG","ANGPTL","ANGPT","CCL","IL16","KLK","SPP1","PLAU","IFNG","CSF","PROS")#top20
  ),
  En_EMTEn = list(
    "COX4I2" = c("EGF","PDGF","CCL","CSF")#all
  ),
  Ep_imEp = list(
    "CAV1" = c("EGF","VEGF","CXCL","IGF","IGFBP","FASLG","SPP1","CSF"),#top20
    "DUOX2" = c("FGF"),#all
    "GDF15" = c("ANGPTL","BMP","CXCL","FGF","VEGF","GAS","GDF","IGF","IL6","NGF","PDGF","SEMA3","EGF","FASLG","ANGPT","CCL","IL16","TNF","IFNG","BAG","SPP1","IGFBP","CSF","EDN"),#all
    "LCN2" = c("GRN","CXCL","IL6","IL16","IGF","IGFBP","FASLG","OSM","CCL","TNF","VEGF","CSF","IFNG","MIF")#all
  ),
  CAF_nCAF = list(
    "NDRG1" = c("IGFBP","KLK","ANGPT","CCL","FGF","PDGF","VEGF","FASLG","TNF","BAG","BMP","CXCL","CSF","IFNG","EDN","PLAU"),#all
    "CDO1" = c("IGF","EGF","PDGF"),#all
    "EGFR" = c("ANGPTL","BMP","CXCL","FGF","VEGF","GAS","GDF","IGF","IL6","NGF","PDGF","SEMA3","EGF","FASLG","TNF","CCL","IL16","ANGPT","EDN","CSF"),#top20
    "TXN" = c("EGF","PDGF","CCL","IGFBP","IL6","CXCL","GRN","FGF","TNF","FASLG","PROS","SPP1","CSF","KLK","IFNG","PLAU"),#top20
    "RARRES2" = c("EGF","VEGF","IGF","PDGF","PROS"),#all
    "TNFAIP3" = c("TNF","GRN","CXCL","IL6","SPP1","CCL","FASLG","IFNG","PDGF","EGF","ANGPTL","BMP","FGF","VEGF","GAS","GDF","IGF","SEMA3","NGF","BAG","KLK"),#top20
    "DPEP1" = c("IL6","EGF","BAG"),#all
    "NOX4" = c("TNF","IL6","EGF","VEGF","CCL","ANGPTL","BMP","CXCL","FGF","GAS","GDF","IGF","NGF","PDGF","SEMA3","FASLG","ANGPT","CSF","KLK","MIF","IFNG")#all
  ),
  PSC_iPSC = list(
    "IL6" = c("ANGPTL","BMP","CXCL","FGF","VEGF","GAS","GDF","IGF","IL6","NGF","PDGF","SEMA3","FASLG","TNF","CCL","ANGPT","IL16","EGF","CSF","SPP1","IFNG","OSM","MIF"),#top20
    "CDKN1A" = c("EGF","FGF","VEGF","IGF","PDGF","MIF","CCL","IL6","CSF","EDN","OSM"),#top20
    "KLF2" = c("EGF","IL6","VEGF","FGF","IGF","CXCL"),#all
    "GABARAPL1" = c("CXCL","CCL","PLAU","IL6","IGF","IGFBP","IFNG"),#all
    "FABP4" = c("ANGPTL","TNF","IGF","IGFBP","BMP","PDGF","FGF","SPP1","IL6","EGF","BAG","CCL","CXCL","PLAU","IFNG","VEGF","NGF","CSF","EDN")#all
  )
)

all_autocrine_networks <- list()

#Function to build autocrine regulatory network
build_autocrine_network <- function(cell_type, actual_col_name, diff_genes) {
  g <- graph.empty(directed = TRUE)
  # step1. Add cell node
  cell_node_name <- cell_type
  g <- add_vertices(g, 1, name = cell_node_name, type = "cell")
  # step2. Add gene nodes
  gene_nodes <- paste0("gene:", make.names(diff_genes))
  g <- add_vertices(g, length(gene_nodes), name = gene_nodes, type = "gene")
  # step3. Add cell-gene edges
  expr_data <- AverageExpression(
    combined_ferro, 
    assays = "RNA", 
    group.by = "full_annotation",
    features = diff_genes
  )$RNA
  cell_expr <- expr_data[, actual_col_name, drop = FALSE]
  
  for(gene in diff_genes) {
    safe_gene <- make.names(gene)
    if(gene %in% rownames(cell_expr)) {
      expr_value <- cell_expr[gene, 1]
      if(expr_value > 0) {
        g <- add_edges(g, 
                       c(cell_node_name, paste0("gene:", safe_gene)),
                       weight = expr_value,
                       type = "expression")
      }
    }
  }
  # step4. Add signaling pathway nodes
  test_pathways <- cellchat_ferro@netP$pathways
  pathway_nodes <- paste0("pathway:", test_pathways)
  g <- add_vertices(g, length(pathway_nodes), name = pathway_nodes, type = "pathway")
  # step5. Add gene-signaling edges
  mapping_key <- gsub(" ", "_", cell_type)
  if(mapping_key %in% names(gene_pathway_mappings)) {
    current_mapping <- gene_pathway_mappings[[mapping_key]]
    
    for(gene in names(current_mapping)) {
      safe_gene <- make.names(gene)
      gene_node <- paste0("gene:", safe_gene)
      if(gene_node %in% V(g)$name) {
        pathways <- current_mapping[[gene]]
        for(pathway in pathways) {
          pathway_node <- paste0("pathway:", pathway)
          if(pathway_node %in% V(g)$name) {
            g <- add_edges(g, 
                           c(gene_node, pathway_node),
                           weight = 1,
                           type = "participation")
          }
        }
      }
    }
  }
  # step6. Add autoregulatory edges (signaling pathways back to the same cell)
  for(pathway in test_pathways) {
    prob_mat <- cellchat_ferro@netP$prob[,, pathway]
    
    if(cell_type %in% rownames(prob_mat) && cell_type %in% colnames(prob_mat)) {
      self_prob <- prob_mat[cell_type, cell_type]
      if(!is.na(self_prob) && self_prob > 0) {
        g <- add_edges(g, 
                       c(paste0("pathway:", pathway), cell_node_name),
                       weight = self_prob,
                       type = "communication")
      }
    }
  }
  return(g)
}

#Build autocrine regulatory networks for all cell types
for(i in seq_along(target_cell_types)) {
  cell_type <- target_cell_types[i]
  actual_col_name <- actual_col_names[i]
  diff_genes <- diff_genes_list[[i]]
  cat("Building autocrine network for:", cell_type, "\n")
  g <- build_autocrine_network(cell_type, actual_col_name, diff_genes)
  
  #Remove isolated signaling pathway nodes
  isolated_pathways <- which(
    V(g)$type == "pathway" & 
      degree(g, v = V(g), mode = "all") == 0
  )
  g <- delete.vertices(g, isolated_pathways)
  all_autocrine_networks[[cell_type]] <- g
  
  #Visualize complete autocrine regulatory network
  V(g)$label <- gsub("^gene:|^pathway:", "", V(g)$name)
  vertex_colors <- c(
    "cell" = "#66C2A5",        
    "pathway" = "#FC8D62",     
    "gene" = "#8DA0CB"        
  )
  vertex_sizes <- c(
    "cell" = 20,
    "pathway" = 14,
    "gene" = 8
  )
  output_file <- paste0("Autocrine_cellgenesignalingnetwork_", 
                        gsub("[ -]", "_", cell_type), ".pdf")
  output_path <- paste0(output_dir, output_file)
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  set.seed(123)
  layout <- layout_with_dh(g)
  pdf(output_path, width = 10, height = 6)
  plot(g,
       layout = layout,
       vertex.size = vertex_sizes[V(g)$type],
       vertex.color = vertex_colors[V(g)$type],
       vertex.label = V(g)$label,
       vertex.label.cex = 0.5,
       vertex.label.color = "black",
       vertex.frame.color = NA,
       edge.arrow.size = 0.4,
       edge.width = 2,  
       edge.color = ifelse(E(g)$type == "expression", "#82cdb5",
                           ifelse(E(g)$type == "participation", "#b5c2dd", "#fda988")),
       main = paste("Autocrine Network for", cell_type),
       main.cex = 0.5)
  legend("bottomright",
         legend = c("Cell", "Gene", "Signaling Pathway", 
                    "Expression", "Participation", "Communication"),
         col = c(vertex_colors["cell"], vertex_colors["gene"], vertex_colors["pathway"],
                 "#82cdb5", "#b5c2dd", "#fda988"),
         pch = c(19, 19, 19, NA, NA, NA),
         lty = c(NA, NA, NA, 1, 1, 1),
         pt.cex = c(2, 0.8, 1, NA, NA, NA),
         bty = "n",
         cex = 0.8)
  dev.off()
  assign(paste0("g_", gsub("[ -]", "_", cell_type), "_autocrine"), g)
}


###Extract autocrine closed-loop subnetworks for each cell type
#Function to extract closed loops
extract_cell_gene_pathway_loops <- function(g) {
  cell_nodes <- V(g)[V(g)$type == "cell"]$name
  all_loops <- list()
  
  for(cell in cell_nodes) {
    # step1. Find all genes expressed by this cell (expression relationships)
    expressed_genes <- names(neighbors(g, cell, mode = "out"))
    expressed_genes <- expressed_genes[startsWith(expressed_genes, "gene:")]
    # step2. For each gene, identify the pathways it participates in
    for(gene in expressed_genes) {
      pathways <- names(neighbors(g, gene, mode = "out"))
      pathways <- pathways[startsWith(pathways, "pathway:")]
      # step3. For each pathway, check if there is communication back to the original cell
      for(pathway in pathways) {
        receiving_cells <- names(neighbors(g, pathway, mode = "out"))
        receiving_cells <- receiving_cells[receiving_cells == cell] 
        
        if(length(receiving_cells) > 0) {
          # step4. Identify complete closed loop
          loop <- list(
            cell_gene = c(cell, gene),
            gene_pathway = c(gene, pathway),
            pathway_cell = c(pathway, cell)
          )
          all_loops[[length(all_loops)+1]] <- loop
        }
      }
    }
  }
  
  if(length(all_loops) == 0) {
    return(list(graph = graph.empty(directed = TRUE), loops = list()))
  }
  # Collect all nodes and edges involved in closed loops
  loop_nodes <- unique(unlist(lapply(all_loops, function(x) unlist(x))))
  loop_edges <- do.call(rbind, lapply(all_loops, function(x) {
    rbind(x$cell_gene, x$gene_pathway, x$pathway_cell)
  }))
  # Create subgraph
  sub_g <- induced_subgraph(g, V(g)[V(g)$name %in% loop_nodes])
  # Retain only edges that are part of closed loops
  edges_to_keep <- apply(loop_edges, 1, function(x) {
    get.edge.ids(sub_g, x)
  })
  edges_to_keep <- edges_to_keep[edges_to_keep != 0] 
  sub_g <- subgraph.edges(sub_g, edges_to_keep)
  return(list(graph = sub_g, loops = all_loops))
}

#Function to build autocrine network with closed loops (only retaining components that form closed loops)
build_autocrine_closed_loop_network <- function(cell_type, actual_col_name, diff_genes) {
  g <- graph.empty(directed = TRUE)
  # step1. Add cell node
  cell_node_name <- cell_type
  g <- add_vertices(g, 1, name = cell_node_name, type = "cell")
  # step2. Add gene nodes
  gene_nodes <- paste0("gene:", make.names(diff_genes))
  g <- add_vertices(g, length(gene_nodes), name = gene_nodes, type = "gene")
  # step3. Add cell-gene edges
  expr_data <- AverageExpression(
    combined_ferro, 
    assays = "RNA", 
    group.by = "full_annotation",
    features = diff_genes
  )$RNA
  cell_expr <- expr_data[, actual_col_name, drop = FALSE]
  
  for(gene in diff_genes) {
    safe_gene <- make.names(gene)
    if(gene %in% rownames(cell_expr)) {
      expr_value <- cell_expr[gene, 1]
      if(expr_value > 0) {
        g <- add_edges(g, 
                       c(cell_node_name, paste0("gene:", safe_gene)),
                       weight = expr_value,
                       type = "expression")
      }
    }
  }
  # step4. Add signaling pathway nodes
  test_pathways <- cellchat_ferro@netP$pathways
  pathway_nodes <- paste0("pathway:", test_pathways)
  g <- add_vertices(g, length(pathway_nodes), name = pathway_nodes, type = "pathway")
  # step5. Add gene-signaling edges (only retain those that can form closed loops)
  mapping_key <- gsub(" ", "_", cell_type)
  if(mapping_key %in% names(gene_pathway_mappings)) {
    current_mapping <- gene_pathway_mappings[[mapping_key]]
    
    for(gene in names(current_mapping)) {
      safe_gene <- make.names(gene)
      gene_node <- paste0("gene:", safe_gene)
      if(gene_node %in% V(g)$name) {
        pathways <- current_mapping[[gene]]
        for(pathway in pathways) {
          pathway_node <- paste0("pathway:", pathway)
          if(pathway_node %in% V(g)$name) {
            # Check if this pathway can form an autocrine closed loop
            prob_mat <- cellchat_ferro@netP$prob[,, pathway]
            if(cell_type %in% rownames(prob_mat) && cell_type %in% colnames(prob_mat)) {
              self_prob <- prob_mat[cell_type, cell_type]
              if(!is.na(self_prob) && self_prob > 0) {
                g <- add_edges(g, 
                               c(gene_node, pathway_node),
                               weight = 1,
                               type = "participation")
              }
            }
          }
        }
      }
    }
  }
  # step6. Add autoregulatory edges (only retain those with upstream gene connections)
  for(pathway in test_pathways) {
    pathway_node <- paste0("pathway:", pathway)
    # Check if this pathway has upstream gene connections
    if(pathway_node %in% V(g)$name && any(incident(g, pathway_node, mode = "in"))) {
      prob_mat <- cellchat_ferro@netP$prob[,, pathway]
      
      if(cell_type %in% rownames(prob_mat) && cell_type %in% colnames(prob_mat)) {
        self_prob <- prob_mat[cell_type, cell_type]
        if(!is.na(self_prob) && self_prob > 0) {
          g <- add_edges(g, 
                         c(pathway_node, cell_node_name),
                         weight = self_prob,
                         type = "communication")
        }
      }
    }
  }
  # step7. Extract closed-loop subgraph
  loop_result <- extract_cell_gene_pathway_loops(g)
  closed_loop_g <- loop_result$graph
  return(list(full_graph = g, closed_loop_graph = closed_loop_g, loops = loop_result$loops))
}

#Build autocrine closed-loop networks for each cell type
for(i in seq_along(target_cell_types)) {
  cell_type <- target_cell_types[i]
  actual_col_name <- actual_col_names[i]
  diff_genes <- diff_genes_list[[i]]
  cat("Building autocrine closed-loop network for:", cell_type, "\n")
  result <- build_autocrine_closed_loop_network(cell_type, actual_col_name, diff_genes)
  
  closed_loop_g <- result$closed_loop_graph
  detected_loops <- result$loops
  
  # Visualize closed-loop subgraph
  if(vcount(closed_loop_g) > 0) {
    cat("Found", length(detected_loops), "complete cell-gene-pathway-cell loops for", cell_type, ":\n")
    
    for(j in seq_along(detected_loops)) {
      loop <- detected_loops[[j]]
      cat(sprintf("Loop %d: %s -> %s -> %s -> %s\n", 
                  j,
                  loop$cell_gene[1],
                  gsub("^gene:", "", loop$cell_gene[2]),
                  gsub("^pathway:", "", loop$gene_pathway[2]),
                  loop$pathway_cell[2]))
    }
    
    V(closed_loop_g)$label <- gsub("^gene:|^pathway:", "", V(closed_loop_g)$name)
    vertex_colors <- c("cell" = "#66C2A5", "pathway" = "#FC8D62", "gene" = "#8DA0CB")
    vertex_sizes <- c("cell" = 20, "pathway" = 14, "gene" = 8)
    
    output_file <- paste0("Autocrine_closedloopnetwork_", 
                          gsub("[ -]", "_", cell_type), ".pdf")
    output_path <- paste0(output_dir, output_file)
    
    set.seed(123)
    layout_sub <- layout_with_dh(closed_loop_g)
    pdf(output_path, width = 10, height = 6)
    plot(closed_loop_g,
         layout = layout_sub,
         vertex.size = vertex_sizes[V(closed_loop_g)$type],
         vertex.color = vertex_colors[V(closed_loop_g)$type],
         vertex.label = V(closed_loop_g)$label,
         vertex.label.cex = 0.5,
         vertex.label.color = "black",
         vertex.frame.color = NA,
         edge.arrow.size = 0.4,
         edge.width = 2,  
         edge.color = ifelse(E(closed_loop_g)$type == "expression", "#82cdb5",
                             ifelse(E(closed_loop_g)$type == "participation", "#b5c2dd", "#fda988")),
         main = paste("Autocrine Closed-Loop Network for", cell_type),
         main.cex = 0.5)
    legend("bottomright",
           legend = c("Cell", "Gene", "Signaling Pathway", 
                      "Expression", "Participation", "Communication"),
           col = c(vertex_colors["cell"], vertex_colors["gene"], vertex_colors["pathway"],
                   "#82cdb5", "#b5c2dd", "#fda988"),
           pch = c(19, 19, 19, NA, NA, NA),
           lty = c(NA, NA, NA, 1, 1, 1),
           pt.cex = c(2, 0.8, 1, NA, NA, NA),
           bty = "n",
           cex = 0.8)
    dev.off()
    
    assign(paste0("g_", gsub("[ -]", "_", cell_type), "_autocrine_closed"), closed_loop_g)
  } else {
    message("No closed loops found in the network for ", cell_type, ".")
  }
}



###Paracrine regulatory networks for each cell type
all_networks <- list()

#Function to build Paracrine regulatory network
build_paracrine_network <- function(cell_type, actual_col_name, diff_genes) {
  g <- graph.empty(directed = TRUE)
  #step1. Add cell node
  cell_node_name <- cell_type  
  g <- add_vertices(g, 1, name = cell_node_name, type = "cell")
  #step2. Add gene nodes
  gene_nodes <- paste0("gene:", make.names(diff_genes))  
  g <- add_vertices(g, length(gene_nodes), name = gene_nodes, type = "gene")
  #step3. Add cell-gene edges
  expr_data <- AverageExpression(
    combined_ferro, 
    assays = "RNA", 
    group.by = "full_annotation",
    features = diff_genes
  )$RNA
  cell_expr <- expr_data[, actual_col_name, drop = FALSE]
  for(gene in diff_genes) {
    safe_gene <- make.names(gene)  
    if(gene %in% rownames(cell_expr)) {
      expr_value <- cell_expr[gene, 1]
      if(expr_value > 0) {
        g <- add_edges(g, 
                       c(cell_node_name, paste0("gene:", safe_gene)),
                       weight = expr_value,
                       type = "expression")
      }
    }
  }
  #step4. Add signaling pathway nodes
  test_pathways <- cellchat_ferro@netP$pathways
  pathway_nodes <- paste0("pathway:", test_pathways)
  g <- add_vertices(g, length(pathway_nodes), name = pathway_nodes, type = "pathway")
  #step5. Add gene-signaling edges
  mapping_key <- gsub(" ", "_", cell_type)  
  if(mapping_key %in% names(gene_pathway_mappings)) {
    current_mapping <- gene_pathway_mappings[[mapping_key]]
    
    for(gene in names(current_mapping)) {
      safe_gene <- make.names(gene)  
      gene_node <- paste0("gene:", safe_gene)
      if(gene_node %in% V(g)$name) {
        pathways <- current_mapping[[gene]]
        for(pathway in pathways) {
          pathway_node <- paste0("pathway:", pathway)
          if(pathway_node %in% V(g)$name) {
            g <- add_edges(g, 
                           c(gene_node, pathway_node),
                           weight = 1,
                           type = "participation")
          }
        }
      }
    }
  }
  #step6. Add paracrine regulatory edges
  for(pathway in test_pathways) {
    prob_mat <- cellchat_ferro@netP$prob[,, pathway]
    
    if(cell_type %in% rownames(prob_mat)) {
      target_cells <- colnames(prob_mat)[prob_mat[cell_type, ] > 0 & 
                                           colnames(prob_mat) != cell_type]
      
      for(target_cell in target_cells) {
        prob_value <- prob_mat[cell_type, target_cell]
        if(!is.na(prob_value) && prob_value > 0) {
          # Ensure target cell node exists
          if(!target_cell %in% V(g)$name) {
            g <- add_vertices(g, 1, name = target_cell, type = "cell")
          }
          # Add edge
          g <- add_edges(g, 
                         c(paste0("pathway:", pathway), target_cell),
                         weight = prob_value,
                         type = "communication")
        }
      }
    }
  }
  return(g)
}

#Build paracrine regulatory networks for all cell types
for(i in seq_along(target_cell_types)) {
  cell_type <- target_cell_types[i]
  actual_col_name <- actual_col_names[i]
  diff_genes <- diff_genes_list[[i]]
  
  cat("Building paracrine network for:", cell_type, "\n")
  g <- build_paracrine_network(cell_type, actual_col_name, diff_genes)
  
  #Remove isolated signaling pathway nodes
  isolated_pathways <- which(
    V(g)$type == "pathway" & 
      degree(g, v = V(g), mode = "all") == 0
  )
  g <- delete.vertices(g, isolated_pathways)
  
  all_networks[[cell_type]] <- g
  
  #Visualization
  V(g)$label <- gsub("^gene:|^pathway:", "", V(g)$name)
  vertex_colors <- c(
    "cell" = "#66C2A5",        
    "pathway" = "#FC8D62",     
    "gene" = "#8DA0CB"        
  )
  vertex_sizes <- c(
    "cell" = 20,
    "pathway" = 14,
    "gene" = 8
  )
  
  output_file <- paste0("Paracrine_cellgenesignalingnetwork_", 
                        gsub("[ -]", "_", cell_type), ".pdf")
  output_path <- paste0(output_dir, output_file)
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  set.seed(123)
  layout <- layout_with_dh(g)
  pdf(output_path, width = 10, height = 6)
  plot(g,
       layout = layout,
       vertex.size = vertex_sizes[V(g)$type],
       vertex.color = vertex_colors[V(g)$type],
       vertex.label = V(g)$label,
       vertex.label.cex = 0.5,
       vertex.label.color = "black",
       vertex.frame.color = NA,
       edge.arrow.size = 0.4,
       edge.width = 2,  
       edge.color = ifelse(E(g)$type == "expression", "#82cdb5",
                           ifelse(E(g)$type == "participation", "#b5c2dd", "#fda988")),
       main = paste("Paracrine Network for", cell_type),
       main.cex = 0.5)
  legend("bottomright",
         legend = c("Cell", "Gene", "Signaling Pathway", 
                    "Expression", "Participation", "Communication"),
         col = c(vertex_colors["cell"], vertex_colors["gene"], vertex_colors["pathway"],
                 "#82cdb5", "#b5c2dd", "#fda988"),
         pch = c(19, 19, 19, NA, NA, NA),
         lty = c(NA, NA, NA, 1, 1, 1),
         pt.cex = c(2, 0.8, 1, NA, NA, NA),
         bty = "n",
         cex = 0.8)
  dev.off()
  assign(paste0("g_", gsub("[ -]", "_", cell_type)), g)
}


###Extract Paracrine Closed-Loop Subnetworks for Each Cell Type
all_networks <- list()
#Function to build paracrine network with closed loops (only retaining components that form closed loops)
build_paracrine_network <- function(cell_type, actual_col_name, diff_genes) {
  g <- graph.empty(directed = TRUE)
  #step1. Add cell node
  cell_node_name <- cell_type
  g <- add_vertices(g, 1, name = cell_node_name, type = "cell")
  #step2. Add gene nodes
  gene_nodes <- paste0("gene:", make.names(diff_genes))
  g <- add_vertices(g, length(gene_nodes), name = gene_nodes, type = "gene")
  #step3. Add cell-gene edges
  expr_data <- AverageExpression(
    combined_ferro, 
    assays = "RNA", 
    group.by = "full_annotation",
    features = diff_genes
  )$RNA
  cell_expr <- expr_data[, actual_col_name, drop = FALSE]
  
  for(gene in diff_genes) {
    safe_gene <- make.names(gene)
    if(gene %in% rownames(cell_expr)) {
      expr_value <- cell_expr[gene, 1]
      if(expr_value > 0) {
        g <- add_edges(g, 
                       c(cell_node_name, paste0("gene:", safe_gene)),
                       weight = expr_value,
                       type = "expression")
      }
    }
  }
  #step4. Add signaling pathway nodes
  test_pathways <- cellchat_ferro@netP$pathways
  pathway_nodes <- paste0("pathway:", test_pathways)
  g <- add_vertices(g, length(pathway_nodes), name = pathway_nodes, type = "pathway")
  
  #step5. Add gene-signaling edges (only retain those that can connect to other cells)
  mapping_key <- gsub(" ", "_", cell_type)
  if(mapping_key %in% names(gene_pathway_mappings)) {
    current_mapping <- gene_pathway_mappings[[mapping_key]]
    
    for(gene in names(current_mapping)) {
      safe_gene <- make.names(gene)
      gene_node <- paste0("gene:", safe_gene)
      if(gene_node %in% V(g)$name) {
        pathways <- current_mapping[[gene]]
        for(pathway in pathways) {
          pathway_node <- paste0("pathway:", pathway)
          if(pathway_node %in% V(g)$name) {
            # Check if this pathway can connect to other cells
            prob_mat <- cellchat_ferro@netP$prob[,, pathway]
            if(cell_type %in% rownames(prob_mat)) {
              target_cells <- colnames(prob_mat)[prob_mat[cell_type, ] > 0 & 
                                                   colnames(prob_mat) != cell_type]
              if(length(target_cells) > 0) {
                g <- add_edges(g, 
                               c(gene_node, pathway_node),
                               weight = 1,
                               type = "participation")
              }
            }
          }
        }
      }
    }
  }
  #step6. Add paracrine regulatory edges (only retain those with upstream gene connections)
  for(pathway in test_pathways) {
    pathway_node <- paste0("pathway:", pathway)
    # Check if this pathway has upstream gene connections
    if(pathway_node %in% V(g)$name && any(incident(g, pathway_node, mode = "in"))) {
      prob_mat <- cellchat_ferro@netP$prob[,, pathway]
      
      if(cell_type %in% rownames(prob_mat)) {
        target_cells <- colnames(prob_mat)[prob_mat[cell_type, ] > 0 & 
                                             colnames(prob_mat) != cell_type]
        
        for(target_cell in target_cells) {
          prob_value <- prob_mat[cell_type, target_cell]
          if(!is.na(prob_value) && prob_value > 0) {
            # Ensure target cell node exists
            if(!target_cell %in% V(g)$name) {
              g <- add_vertices(g, 1, name = target_cell, type = "cell")
            }
            # Add edge
            g <- add_edges(g, 
                           c(pathway_node, target_cell),
                           weight = prob_value,
                           type = "communication")
          }
        }
      }
    }
  }
  #step7. Remove isolated nodes (nodes that don't form complete paths)
  #First identify all nodes that are part of complete paths
  valid_nodes <- c()
  for(pathway in test_pathways) {
    pathway_node <- paste0("pathway:", pathway)
    if(pathway_node %in% V(g)$name) {
      #Check if there are upstream genes and downstream cells
      in_edges <- incident(g, pathway_node, mode = "in")
      out_edges <- incident(g, pathway_node, mode = "out")
      if(length(in_edges) > 0 && length(out_edges) > 0) {
        #Add upstream genes and cells
        gene_nodes <- ends(g, in_edges)[,1]
        valid_nodes <- c(valid_nodes, gene_nodes, pathway_node)
        #Add downstream cells
        cell_nodes <- ends(g, out_edges)[,2]
        valid_nodes <- c(valid_nodes, cell_nodes)
      }
    }
  }
  #Keep all valid nodes and their connected cell nodes
  valid_nodes <- unique(c(cell_node_name, valid_nodes))
  nodes_to_remove <- setdiff(V(g)$name, valid_nodes)
  g <- delete.vertices(g, nodes_to_remove)
  
  return(g)
}


#Build paracrine closed-loop networks for all cell types
for(i in seq_along(target_cell_types)) {
  cell_type <- target_cell_types[i]
  actual_col_name <- actual_col_names[i]
  diff_genes <- diff_genes_list[[i]]
  
  cat("Building paracrine network for:", cell_type, "\n")
  g <- build_paracrine_network(cell_type, actual_col_name, diff_genes)
  
  #Remove isolated signaling pathway nodes
  isolated_pathways <- which(
    V(g)$type == "pathway" & 
      degree(g, v = V(g), mode = "all") == 0
  )
  g <- delete.vertices(g, isolated_pathways)
  
  all_networks[[cell_type]] <- g
  
  #Visualize paracrine closed-loop subgraph
  V(g)$label <- gsub("^gene:|^pathway:", "", V(g)$name)
  vertex_colors <- c(
    "cell" = "#66C2A5",        
    "pathway" = "#FC8D62",     
    "gene" = "#8DA0CB"        
  )
  vertex_sizes <- c(
    "cell" = 20,
    "pathway" = 14,
    "gene" = 8
  )
  
  output_file <- paste0("Paracrine_closedloopnetwork_", 
                        gsub("[ -]", "_", cell_type), ".pdf")
  output_path <- paste0(output_dir, output_file)
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  set.seed(123)
  layout <- layout_with_dh(g)
  pdf(output_path, width = 10, height = 6)
  plot(g,
       layout = layout,
       vertex.size = vertex_sizes[V(g)$type],
       vertex.color = vertex_colors[V(g)$type],
       vertex.label = V(g)$label,
       vertex.label.cex = 0.5,
       vertex.label.color = "black",
       vertex.frame.color = NA,
       edge.arrow.size = 0.4,
       edge.width = 2,  
       edge.color = ifelse(E(g)$type == "expression", "#82cdb5",
                           ifelse(E(g)$type == "participation", "#b5c2dd", "#fda988")),
       main = paste("Paracrine Network for", cell_type),
       main.cex = 0.5)
  legend("bottomright",
         legend = c("Cell", "Gene", "Signaling Pathway", 
                    "Expression", "Participation", "Communication"),
         col = c(vertex_colors["cell"], vertex_colors["gene"], vertex_colors["pathway"],
                 "#82cdb5", "#b5c2dd", "#fda988"),
         pch = c(19, 19, 19, NA, NA, NA),
         lty = c(NA, NA, NA, 1, 1, 1),
         pt.cex = c(2, 0.8, 1, NA, NA, NA),
         bty = "n",
         cex = 0.8)
  dev.off()
  
  assign(paste0("g_", gsub("[ -]", "_", cell_type)), g)
}


