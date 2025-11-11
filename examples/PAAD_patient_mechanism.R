######Ferroptosis mechanism in PAAD patients
######################################################



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("tidyverse","pheatmap","clusterProfiler","org.Hs.eg.db",
                       "ggrepel","tidytext","ggplot2","dplyr","patchwork",
                       "readr","tibble","ggnewscale","msigdbr","enrichplot",
                       "rmda","DESeq2","mgcv","rstatix","PMCMRplus","vcd",
                       "ggpubr","coin","reshape2"))

library(tidyverse)    
library(pheatmap)     
library(clusterProfiler)  
library(org.Hs.eg.db) 
library(ggrepel)      
library(tidytext)  
library(ggplot2)
library(dplyr)
library(patchwork)
library(readr)
library(tibble)
library(ggnewscale)
library(msigdbr)
library(enrichplot)
library(rmda)
library(DESeq2)
library(mgcv)
library(rstatix)
library(PMCMRplus)
library(vcd)       
library(ggpubr)
library(coin)
library(reshape2)


output_dir <- "my_results/mechanism"



#Ferroptosis characteristics of patients across different tumor stages
#DEGs in patients across different tumor stages
stage1_samples <- read_csv("my_data/Stage1_id.csv")  
stage1_deg_files <- list.files("my_data/Stage1_DEG", pattern = "\\.csv$", full.names = TRUE)
stage1_deg_all <- map_dfr(stage1_deg_files, ~ {
  read_csv(.x) %>%
    mutate(patient_id = str_extract(basename(.x), "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")) 
})
stage1_deg_sig <- stage1_deg_all %>%
  filter(change != "Not Change")
write.csv(stage1_deg_sig, file = file.path(output_dir, "deg_sig_Stage1.csv"), quote = F)

stage2_samples <- read_csv("my_data/Stage2_id.csv")  
stage2_deg_files <- list.files("my_data/Stage2_DEG", pattern = "\\.csv$", full.names = TRUE)
stage2_deg_all <- map_dfr(stage2_deg_files, ~ {
  read_csv(.x) %>%
    mutate(patient_id = str_extract(basename(.x), "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")) 
})
stage2_deg_sig <- stage2_deg_all %>%
  filter(change != "Not Change")
write.csv(stage2_deg_sig, file = file.path(output_dir, "deg_sig_Stage2.csv"), quote = F)

stage3_samples <- read_csv("my_data/Stage3_id.csv")  
stage3_deg_files <- list.files("my_data/Stage3_DEG", pattern = "\\.csv$", full.names = TRUE)
stage3_deg_all <- map_dfr(stage3_deg_files, ~ {
  read_csv(.x) %>%
    mutate(patient_id = str_extract(basename(.x), "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")) 
})
stage3_deg_sig <- stage3_deg_all %>%
  filter(change != "Not Change")
write.csv(stage3_deg_sig, file = file.path(output_dir, "deg_sig_Stage3.csv"), quote = F)

stage4_samples <- read_csv("my_data/Stage4_id.csv")  
stage4_deg_files <- list.files("my_data/Stage4_DEG", pattern = "\\.csv$", full.names = TRUE)
stage4_deg_all <- map_dfr(stage4_deg_files, ~ {
  read_csv(.x) %>%
    mutate(patient_id = str_extract(basename(.x), "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")) 
})
stage4_deg_sig <- stage4_deg_all %>%
  filter(change != "Not Change")
write.csv(stage4_deg_sig, file = file.path(output_dir, "deg_sig_Stage4.csv"), quote = F)


#High-frequency consensus DEGs in patients across different tumor stages
deg_summary <- bind_rows(
  stage1_deg_sig %>% mutate(group = "Stage1"),
  stage2_deg_sig %>% mutate(group = "Stage2"),
  stage3_deg_sig %>% mutate(group = "Stage3"),
  stage4_deg_sig %>% mutate(group = "Stage4")
) %>%
  count(group, change, name = "count") %>%
  complete(group, change, fill = list(count = 0)) 
write.csv(deg_summary, file = file.path(output_dir, "deg_summary_Stage.csv"), quote = F)

#Extract consensus DEGs that appear frequently in each group (consistently changed in at least 50% of patients)
get_consensus_genes <- function(deg_data, min_patients_ratio = 0.5) {
  deg_data %>%
    group_by(Gene, change) %>%
    summarise(n_patients = n_distinct(patient_id), .groups = "drop") %>%
    filter(n_patients >= max(n_patients) * min_patients_ratio) %>%
    arrange(desc(n_patients))
}
stage1_consensus <- get_consensus_genes(stage1_deg_sig)
stage2_consensus <- get_consensus_genes(stage2_deg_sig)
stage3_consensus <- get_consensus_genes(stage3_deg_sig)
stage4_consensus <- get_consensus_genes(stage4_deg_sig)

#Merge high-frequency consensus genes from all groups
all_consensus <- bind_rows(
  stage1_consensus %>% mutate(group = "Stage1"),
  stage2_consensus %>% mutate(group = "Stage2"),
  stage3_consensus %>% mutate(group = "Stage3"),
  stage4_consensus %>% mutate(group = "Stage4")
) 
write.csv(all_consensus, file = file.path(output_dir, "consensusgene_all_Stage.csv"), quote = F)

#Extract top 20 high-frequency consensus genes from each group
top_genes <- bind_rows(
  stage1_consensus %>% head(20) %>% mutate(group = "Stage1"),
  stage2_consensus %>% head(20) %>% mutate(group = "Stage2"),
  stage3_consensus %>% head(20) %>% mutate(group = "Stage3"),
  stage4_consensus %>% head(20) %>% mutate(group = "Stage4")
) %>% distinct(Gene, change, n_patients, group)
write.csv(top_genes, file = file.path(output_dir, "consensusgene_top_Stage.csv"), quote = F)

#Draw grouped bar chart
pdf(file.path(output_dir, "top_consensusgene_barchat_Stage.pdf"), width = 14, height = 10)
stage_labels <- c(
  "Stage1" = "Stage I",
  "Stage2" = "Stage II",
  "Stage3" = "Stage III",
  "Stage4" = "Stage IV"
)
ggplot(top_genes, aes(x = reorder_within(Gene, n_patients, group), 
                      y = n_patients, fill = change)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +  
  scale_fill_manual(values = c("Up" = "#e12216", "Down" = "#154996"),  
                    name = "Expression Change",  
                    labels = c("Down", "Up")) +  
  scale_x_reordered() +
  facet_wrap(~ group, scales = "free_y", ncol = 2,
             labeller = labeller(group = stage_labels)) +  
  coord_flip() +
  labs(x = "Gene", 
       y = "Number of Patients") +  
  theme_minimal(base_size = 12) +  
  theme(
    strip.text = element_text(face = "bold", size = 15, hjust = 0.5), 
    panel.grid.major.y = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.title = element_text(face = "bold", size = 21),  
    axis.text = element_text(color = "black", size = 16),  
    legend.position = "top",  
    legend.title = element_text(face = "bold"),  
    panel.spacing = unit(1, "lines"),  
    plot.margin = margin(15, 15, 15, 15)  
  )
dev.off()


#Heatmap of top DEGs in patients across different tumor stages
#Load original expression matrix
expr_matrix <- read_csv("my_data/Experiencecurve_PAAD_matrix.csv")
colnames(expr_matrix) <- ifelse(
  str_detect(colnames(expr_matrix), "^TCGA\\."),  
  str_replace_all(colnames(expr_matrix), "\\.", "-"),  
  colnames(expr_matrix)  
)

valid_matrix<-read_csv("my_data/ferroptosis_scores_PAAD.csv")
colnames(valid_matrix) <- ifelse(
  str_detect(colnames(valid_matrix), "^TCGA\\."),  
  str_replace_all(colnames(valid_matrix), "\\.", "-"),  
  colnames(valid_matrix)  
)

#Filter original expression matrix for high-frequency consensus genes
expr_matrix <- as_tibble(expr_matrix)
heatmap_data <- expr_matrix %>%
  dplyr::select(Gene, any_of(colnames(valid_matrix))) %>%
  filter(Gene %in% top_genes$Gene) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

annotation_col <- data.frame(
  Group = ifelse(colnames(heatmap_data) %in% stage1_samples$patient_id, "Stage1",
                 ifelse(colnames(heatmap_data) %in% stage2_samples$patient_id, "Stage2",
                        ifelse(colnames(heatmap_data) %in% stage3_samples$patient_id, "Stage3", "Stage4"))))
rownames(annotation_col) <- colnames(heatmap_data)

#Draw heatmap
pdf(file.path(output_dir, "top_consensusgene_heatmap_Stage.pdf"), width = 14, height = 10)
annot_col_modified <- annotation_col
annot_col_modified$Group <- factor(annot_col_modified$Group,
                                   levels = c("Stage1", "Stage2", "Stage3", "Stage4"),
                                   labels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
annotation_colors <- list(
  Group = c("Stage1" = "#FFCCCC", "Stage2" = "#FFE6CC", "Stage3" = "#FFFFCC", "Stage4" = "#CCFFCC"))
annot_colors_modified <- annotation_colors
names(annot_colors_modified$Group) <- c("Stage I", "Stage II", "Stage III", "Stage IV")
heatmap_colors <- colorRampPalette(c("#154996", "#f7f7f7", "#e12216"))(100)
pheatmap(
  heatmap_data,
  scale = "row",
  color = heatmap_colors,
  annotation_col = annot_col_modified,
  annotation_colors = annot_colors_modified, 
  show_colnames = FALSE,
  show_rownames = TRUE, 
  fontsize_row = 12, 
  fontsize_col = 12, 
  fontsize = 14, 
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean", 
  clustering_distance_cols = "euclidean",
  clustering_method = "complete", 
  treeheight_row = 60, 
  treeheight_col = 60, 
  legend = TRUE
)
dev.off()



#GO enrichment analysis for patients across different tumor stages
stage1_go_genes <- stage1_deg_sig %>% pull(Gene)
stage1_go_ids <- bitr(stage1_go_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
stage1_go_enrich <- enrichGO(
  gene = stage1_go_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
stage1_transGO <- setReadable(stage1_go_enrich, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(stage1_transGO, file = file.path(output_dir, "GO_Stage1.csv"), quote = F)

stage2_go_genes <- stage2_deg_sig %>% pull(Gene)
stage2_go_ids <- bitr(stage2_go_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
stage2_go_enrich <- enrichGO(
  gene = stage2_go_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
stage2_transGO <- setReadable(stage2_go_enrich, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(stage2_transGO, file = file.path(output_dir, "GO_Stage2.csv"), quote = F)

stage3_go_genes <- stage3_deg_sig %>% pull(Gene)
stage3_go_ids <- bitr(stage3_go_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
stage3_go_enrich <- enrichGO(
  gene = stage3_go_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
stage3_transGO <- setReadable(stage3_go_enrich, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(stage3_transGO, file = file.path(output_dir, "GO_Stage3.csv"), quote = F)

stage4_go_genes <- stage4_deg_sig %>% pull(Gene)
stage4_go_ids <- bitr(stage4_go_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
stage4_go_enrich <- enrichGO(
  gene = stage4_go_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
stage4_transGO <- setReadable(stage4_go_enrich, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(stage4_transGO, file = file.path(output_dir, "GO_Stage4.csv"), quote = F)

#Draw GO plots
p1 <- barplot(stage1_transGO, drop = TRUE, showCategory = 15,label_format = 100) + 
  ggtitle("Stage I") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 21),  
        axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        legend.text = element_text(size = 16),  
        legend.title = element_text(size = 18), 
        panel.spacing = unit(1, "lines"),  
        plot.margin = margin(15, 15, 15, 15))

p2 <- barplot(stage2_transGO, drop = TRUE, showCategory = 15,label_format = 100) + 
  ggtitle("Stage II") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 21),  
        axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        legend.text = element_text(size = 16),  
        legend.title = element_text(size = 18), 
        panel.spacing = unit(1, "lines"),  
        plot.margin = margin(15, 15, 15, 15))

p3 <- barplot(stage3_transGO, drop = TRUE, showCategory = 15,label_format = function(x) str_wrap(x, width = 60)) + 
  ggtitle("Stage III") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 21),  
        axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        legend.text = element_text(size = 16),  
        legend.title = element_text(size = 18), 
        panel.spacing = unit(1, "lines"),  
        plot.margin = margin(15, 15, 15, 15))

p4 <- barplot(stage4_transGO, drop = TRUE, showCategory = 15,label_format = function(x) str_wrap(x, width = 60)) + 
  ggtitle("Stage IV") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 21),  
        axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        legend.text = element_text(size = 16),  
        legend.title = element_text(size = 18), 
        panel.spacing = unit(1, "lines"),  
        plot.margin = margin(15, 15, 15, 15))

combined_plot <- (p1 | p2) / (p3 | p4) +
  plot_annotation(theme = theme(plot.title = element_text(hjust = 0.3, size = 20)))

pdf(file.path(output_dir, "GO_Stage.pdf"), width = 26, height = 16)
print(combined_plot)
dev.off()


#KEGG enrichment analysis for patients across different tumor stages
stage1_kegg <- enrichKEGG(gene = stage1_go_ids$ENTREZID,
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
stage1_kk <- setReadable(stage1_kegg, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(stage1_kk, file = file.path(output_dir, "KEGG_Stage1.csv"), quote = F)

stage2_kegg <- enrichKEGG(gene = stage2_go_ids$ENTREZID,
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
stage2_kk <- setReadable(stage2_kegg, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(stage2_kk, file = file.path(output_dir, "KEGG_Stage2.csv"), quote = F)

stage3_kegg <- enrichKEGG(gene = stage3_go_ids$ENTREZID,
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
stage3_kk <- setReadable(stage3_kegg, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(stage3_kk, file = file.path(output_dir, "KEGG_Stage3.csv"), quote = F)

stage4_kegg <- enrichKEGG(gene = stage4_go_ids$ENTREZID,
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
stage4_kk <- setReadable(stage4_kegg, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(stage4_kk, file = file.path(output_dir, "KEGG_Stage4.csv"), quote = F)

#Draw KEGG plots
p1 <- dotplot(stage1_kegg, showCategory = 15,label_format = 100) + 
  ggtitle("Stage I") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 21),  
        axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        legend.text = element_text(size = 16),  
        legend.title = element_text(size = 18), 
        panel.spacing = unit(1, "lines"),  
        plot.margin = margin(15, 15, 15, 15))

p2 <- dotplot(stage2_kegg, showCategory = 15,label_format = 100) + 
  ggtitle("Stage II") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 21),  
        axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        legend.text = element_text(size = 16),  
        legend.title = element_text(size = 18), 
        panel.spacing = unit(1, "lines"),  
        plot.margin = margin(15, 15, 15, 15))

p3 <- dotplot(stage3_kegg, showCategory = 15,label_format = 100) + 
  ggtitle("Stage III") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 21),  
        axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        legend.text = element_text(size = 16),  
        legend.title = element_text(size = 18), 
        panel.spacing = unit(1, "lines"),  
        plot.margin = margin(15, 15, 15, 15))

p4 <- dotplot(stage4_kegg, showCategory = 15,label_format = 100) + 
  ggtitle("Stage IV") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 21),  
        axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        legend.text = element_text(size = 16),  
        legend.title = element_text(size = 18), 
        panel.spacing = unit(1, "lines"),  
        plot.margin = margin(15, 15, 15, 15))

combined_plot <- (p1 | p2) / (p3 | p4) +
  plot_annotation(theme = theme(plot.title = element_text(hjust = 0.3, size = 16)))

pdf(file.path(output_dir, "KEGG_Stage.pdf"), width = 22, height = 16)
print(combined_plot)
dev.off()





#Ferroptosis characteristics at different levels
#Data preparation
ferro_data <- read_csv("my_data/ferroptosis_scores_PAAD.csv")
clinical_data <- read_csv("my_data/TCGA_PAAD_clinical.csv")

ferro_scores <- ferro_data[4,]
ferro_scores_transposed <- ferro_scores %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "patient_id") %>% 
  mutate(patient_id = gsub("\\.", "-", patient_id)) %>%
  mutate(patient_id = substr(patient_id, 1, 12))
ferro_scores<-ferro_scores_transposed[-1,]
colnames(ferro_scores)[2] <- "Ferroptosis_score"

survival_time <- clinical_data %>% 
  select(patient_id = 1, survival_time = 16)

merged_data <- inner_join(survival_time, ferro_scores, by = "patient_id")

#Grouping based on experience curve inflection points
merged_data$Ferroptosis_score <- as.numeric(merged_data$Ferroptosis_score)
merged_data <- merged_data %>%
  mutate(group = case_when(
    Ferroptosis_score <= 0.5 ~ "Low-Fer",
    Ferroptosis_score >= 0.8 ~ "High-Fer",
    TRUE ~ "Mid-Fer"
  ))
table(merged_data$group)

#Data processing for count and group information
expr_matrix <- read_csv("my_data/Experiencecurve_PAAD_matrix.csv")
colnames(expr_matrix) <- ifelse(
  str_detect(colnames(expr_matrix), "^TCGA\\."),  
  str_replace_all(colnames(expr_matrix), "\\.", "-"),  
  colnames(expr_matrix)  
)
expr_matrix <- as.data.frame(expr_matrix)
rownames(expr_matrix) <- expr_matrix[,1]
expr_matrix <- expr_matrix[,-1]

common_samples <- intersect(colnames(expr_matrix), merged_data$patient_id)
expr_matrix <- expr_matrix[, common_samples]
merged_data <- merged_data %>% filter(patient_id %in% common_samples)
merged_data$group <- factor(merged_data$group, 
                            levels = c("Low-Fer", "Mid-Fer", "High-Fer"))

#DEG analysis
dds <- DESeqDataSetFromMatrix(countData = round(expr_matrix),
                              colData = merged_data,
                              design = ~ group)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

comparisons <- list(
  c("Low-Fer", "Mid-Fer"),
  c("High-Fer", "Mid-Fer"),
  c("High-Fer", "Low-Fer")
)

results_list <- lapply(comparisons, function(comp) {
  res <- results(dds, contrast = c("group", comp[1], comp[2]))
  return(res)
})
names(results_list) <- sapply(comparisons, paste, collapse = "_vs_")

sig_genes <- lapply(results_list, function(res) {
  res_df <- as.data.frame(res) %>% 
    rownames_to_column("gene") %>% 
    filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
    mutate(
      regulation = case_when(
        log2FoldChange > 0 ~ "Up",
        log2FoldChange < 0 ~ "Down",
        TRUE ~ "Not Change"
      )
    ) %>%
    arrange(padj)
  return(res_df)
})

write.csv(bind_rows(sig_genes, .id = "comparison"), file = file.path(output_dir, "differential_genes_results.csv"))

#GO enrichment analysis
all_diff_genes <- bind_rows(sig_genes, .id = "comparison")
go_results <- lapply(sig_genes, function(df) {
  if(nrow(df) > 0) {
    gene_ids <- bitr(df$gene, fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)
    
    enrichGO(gene = gene_ids$ENTREZID,
             OrgDb = org.Hs.eg.db,
             keyType = "ENTREZID",
             ont = "BP",
             pAdjustMethod = "BH",
             qvalueCutoff = 0.05)
  } else {
    return(NULL)
  }
})

pdf(file.path(output_dir, "GO_Low_vs_Mid.pdf"), width = 10, height = 8)
if(!is.null(go_results[["Low-Fer_vs_Mid-Fer"]])) {
  dotplot(go_results[["Low-Fer_vs_Mid-Fer"]], 
          title = "GO Enrichment - Low-Fer vs Mid-Fer",
          showCategory = 15)
}
dev.off()
pdf(file.path(output_dir, "GO_High_vs_Mid.pdf"), width = 10, height = 8)
if(!is.null(go_results[["High-Fer_vs_Mid-Fer"]])) {
  dotplot(go_results[["High-Fer_vs_Mid-Fer"]], 
          title = "GO Enrichment - High-Fer vs Mid-Fer",
          showCategory = 15)
}
dev.off()
pdf(file.path(output_dir, "GO_High_vs_Low.pdf"), width = 10, height = 8)
if(!is.null(go_results[["High-Fer_vs_Low-Fer"]])) {
  dotplot(go_results[["High-Fer_vs_Low-Fer"]], 
          title = "GO Enrichment - High-Fer vs Low-Fer",
          showCategory = 15)
}
dev.off()

#KEGG enrichment analysis
kegg_results <- lapply(sig_genes, function(df) {
  if(nrow(df) > 0) {
    gene_ids <- bitr(df$gene, fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)
    
    enrichKEGG(gene = gene_ids$ENTREZID,
               organism = "hsa",
               pAdjustMethod = "BH",
               qvalueCutoff = 0.05)
  } else {
    return(NULL)
  }
})

pdf(file.path(output_dir, "KEGG_Low_vs_Mid.pdf"), width = 10, height = 8)
if(!is.null(kegg_results[["Low-Fer_vs_Mid-Fer"]])) {
  dotplot(kegg_results[["Low-Fer_vs_Mid-Fer"]], 
          title = "KEGG Enrichment - Low-Fer vs Mid-Fer",
          showCategory = 15)
}
dev.off()
pdf(file.path(output_dir, "KEGG_High_vs_Mid.pdf"), width = 10, height = 8)
if(!is.null(kegg_results[["High-Fer_vs_Mid-Fer"]])) {
  dotplot(kegg_results[["High-Fer_vs_Mid-Fer"]], 
          title = "KEGG Enrichment - High-Fer vs Mid-Fer",
          showCategory = 15)
}
dev.off()
pdf(file.path(output_dir, "KEGG_High_vs_Low.pdf"), width = 10, height = 8)
if(!is.null(kegg_results[["High-Fer_vs_Low-Fer"]])) {
  dotplot(kegg_results[["High-Fer_vs_Low-Fer"]], 
          title = "KEGG Enrichment - High-Fer vs Low-Fer",
          showCategory = 15)
}
dev.off()

#Summary plot
prepare_enrich_data <- function(enrich_result, type, comparison) {
  if(!is.null(enrich_result)) {
    df <- as.data.frame(enrich_result) %>% 
      arrange(p.adjust) %>%  
      head(15) %>%           
      mutate(Type = type,
             Comparison = comparison,
             GeneRatio = sapply(strsplit(GeneRatio, "/"), 
                                function(x) as.numeric(x[1])/as.numeric(x[2])))
    return(df)
  }
  return(NULL)
}
combined_data <- bind_rows(
  prepare_enrich_data(kegg_results[["Low-Fer_vs_Mid-Fer"]], 
                      "KEGG", "KEGG: Low-Fer vs Mid-Fer"),
  prepare_enrich_data(kegg_results[["High-Fer_vs_Mid-Fer"]], 
                      "KEGG", "KEGG: High-Fer vs Mid-Fer"),
  prepare_enrich_data(go_results[["High-Fer_vs_Mid-Fer"]], 
                      "GO", "GO: High-Fer vs Mid-Fer")
) %>%
  mutate(Comparison = factor(Comparison,
                             levels = c("KEGG: Low-Fer vs Mid-Fer",
                                        "KEGG: High-Fer vs Mid-Fer",
                                        "GO: High-Fer vs Mid-Fer")))

create_unified_dotplot <- function(data) {
  min_padj <- min(data$p.adjust, na.rm = TRUE)
  max_count <- max(data$Count, na.rm = TRUE)
  ggplot(data, aes(x = GeneRatio, 
                   y = reorder(Description, GeneRatio))) +
    geom_point(aes(size = Count, color = p.adjust)) +  
    scale_color_gradient(low = "#519D78", high = "#F3FBF2",  
                         limits = c(0, 0.05),  
                         name = "p.adjust") + 
    scale_size_continuous(range = c(1, 10),
                          limits = c(1, max_count),
                          name = "Gene Count") +
    facet_grid(Comparison ~ ., scales = "free_y", space = "free", 
               switch = "y") + 
    labs(x = "Gene Ratio", y = "") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text.y.left = element_blank(),  
      strip.background = element_blank(),  
      axis.text.y = element_text(size = 14, hjust = 1),  
      axis.text.y.right = element_text(size = 14),  
      axis.ticks.y.right = element_line(), 
      axis.line.y.right = element_line(),  
      legend.position.inside = c(0.95, 0.1),  
      legend.justification = c(1, 0),         
      legend.box = "vertical",               
      legend.box.just = "right",             
      legend.direction = "vertical",         
      legend.spacing.y = unit(0.3, "cm"),    
      legend.key.size = unit(0.5, "cm"),     
      legend.text = element_text(size = 10), 
      legend.title = element_text(size = 10) 
    ) +
    scale_y_discrete(position = "right")  
}
final_plot <- create_unified_dotplot(combined_data)
pdf(file.path(output_dir, "Enrichment_GO+KEGG.pdf"), width = 14, height = 12)
print(final_plot)
dev.off()

#GSEA analysis
gsea_data_hl <- results_list[["High-Fer_vs_Low-Fer"]]
gene_rank_hl <- gsea_data_hl$log2FoldChange
names(gene_rank_hl) <- rownames(gsea_data_hl)
gene_rank_hl <- sort(gene_rank_hl, decreasing = TRUE)
ferro_sets <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  filter(grepl("GLUTATHIONE|FATTY_ACID|FERR|LIPID", gs_name, ignore.case = TRUE))
gsea_res_hl <- GSEA(gene_rank_hl, 
                    TERM2GENE = ferro_sets[, c("gs_name", "gene_symbol")],
                    pvalueCutoff = 0.05,
                    eps = 0)
write.csv(as.data.frame(gsea_res_hl), 
          file = file.path(output_dir, "GSEA_High_vs_Low.csv"))
pdf(file.path(output_dir, "GSEA_High_vs_Low.pdf"), width = 12, height = 8)
gsea_colors <- c("#BF242A", "#4D5AAF", "#007B43") 
gseaplot2(gsea_res_hl, 
          geneSetID = 1:min(3, nrow(gsea_res_hl)),
          pvalue_table = TRUE,
          title = "Ferroptosis-related Pathways (High-Fer vs Low-Fer)",
          color = gsea_colors,
          base_size = 18, 
          ES_geom = "line") 
dev.off()

gsea_data_lm <- results_list[["Low-Fer_vs_Mid-Fer"]]
gene_rank_lm <- gsea_data_lm$log2FoldChange
names(gene_rank_lm) <- rownames(gsea_data_lm)
gene_rank_lm <- sort(gene_rank_lm, decreasing = TRUE)
ferro_sets <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  filter(grepl("GLUTATHIONE|FATTY_ACID|FERR|LIPID", gs_name, ignore.case = TRUE))
gsea_res_lm <- GSEA(gene_rank_lm, 
                    TERM2GENE = ferro_sets[, c("gs_name", "gene_symbol")],
                    pvalueCutoff = 0.05,
                    eps = 0)
write.csv(as.data.frame(gsea_res_lm), 
          file = file.path(output_dir, "GSEA_Low_vs_Mid.csv"))
pdf(file.path(output_dir, "GSEA_Low_vs_Mid.pdf"), width = 12, height = 8)
gsea_colors <- c("#BF242A", "#4D5AAF", "#007B43") 
gseaplot2(gsea_res_lm, 
          geneSetID = 1:min(3, nrow(gsea_res_lm)),
          pvalue_table = TRUE,
          title = "Ferroptosis-related Pathways (Low-Fer vs Mid-Fer)",
          color = gsea_colors,
          base_size = 18, 
          ES_geom = "line") 
dev.off()

gsea_data_hm <- results_list[["High-Fer_vs_Mid-Fer"]]
gene_rank_hm <- gsea_data_hm$log2FoldChange
names(gene_rank_hm) <- rownames(gsea_data_hm)
gene_rank_hm <- sort(gene_rank_hm, decreasing = TRUE)
ferro_sets <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  filter(grepl("GLUTATHIONE|FATTY_ACID|FERR|LIPID", gs_name, ignore.case = TRUE))
gsea_res_hm <- GSEA(gene_rank_hm, 
                    TERM2GENE = ferro_sets[, c("gs_name", "gene_symbol")],
                    pvalueCutoff = 0.05,
                    eps = 0)
write.csv(as.data.frame(gsea_res_hm), 
          file = file.path(output_dir, "GSEA_High_vs_Mid.csv"))
pdf(file.path(output_dir, "GSEA_High_vs_Mid.pdf"), width = 12, height = 8)
gsea_colors <- c("#BF242A", "#4D5AAF", "#007B43") 
gseaplot2(gsea_res_hm, 
          geneSetID = 1:min(3, nrow(gsea_res_hm)),
          pvalue_table = TRUE,
          title = "Ferroptosis-related Pathways (High-Fer vs Mid-Fer)",
          color = gsea_colors,
          base_size = 18, 
          ES_geom = "line") 
dev.off()

#Summary plot
gsea_colors <- c("#BF242A", "#4D5AAF", "#007B43")

gsea_data_hl <- results_list[["High-Fer_vs_Low-Fer"]]
gene_rank_hl <- gsea_data_hl$log2FoldChange
names(gene_rank_hl) <- rownames(gsea_data_hl)
gene_rank_hl <- sort(gene_rank_hl, decreasing = TRUE)

gsea_data_hm <- results_list[["High-Fer_vs_Mid-Fer"]]
gene_rank_hm <- gsea_data_hm$log2FoldChange
names(gene_rank_hm) <- rownames(gsea_data_hm)
gene_rank_hm <- sort(gene_rank_hm, decreasing = TRUE)

ferro_sets <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  filter(grepl("GLUTATHIONE|FATTY_ACID|FERR|LIPID", gs_name, ignore.case = TRUE))

gsea_res_hl <- GSEA(gene_rank_hl, 
                    TERM2GENE = ferro_sets[, c("gs_name", "gene_symbol")],
                    pvalueCutoff = 0.05,
                    eps = 0)
gsea_res_hm <- GSEA(gene_rank_hm, 
                    TERM2GENE = ferro_sets[, c("gs_name", "gene_symbol")],
                    pvalueCutoff = 0.05,
                    eps = 0)

pdf(file.path(output_dir, "Enrichment_GSEA.pdf"), width = 12, height = 8)
gseaplot2(gsea_res_hl, 
          geneSetID = 1:min(3, nrow(gsea_res_hl)),
          pvalue_table = FALSE,
          title = "High-Fer vs Low-Fer",
          color = gsea_colors,
          base_size = 16,
          ES_geom = "line",
          subplots = 1)  
gseaplot2(gsea_res_hm, 
          geneSetID = 1:min(3, nrow(gsea_res_hm)),
          pvalue_table = FALSE,
          title = "High-Fer vs Mid-Fer",
          color = gsea_colors,
          base_size = 16,
          ES_geom = "line",
          subplots = 1)
dev.off()





#Decision Curve Analysis
#Calculate ferroptosis grade
ferro_scores<-read_csv("my_data/ferroptosis_scores_PAAD.csv")
ferro_grades_quantile <- ferro_scores %>%
  mutate(
    Ferroptosis_score = as.numeric(Ferroptosis_score),
    
    quantile_cutoffs = list(quantile(Ferroptosis_score, probs = seq(0, 1, 0.1), na.rm = TRUE)),
    
    Ferroptosis_grade = sapply(Ferroptosis_score, function(score) {
      if(is.na(score)) return(1)  
      cutoffs <- unlist(quantile_cutoffs[1])  
      for(i in 10:1) {
        if(score >= cutoffs[i]) return(i)
      }
      return(1)
    })
  ) %>%
  select(-quantile_cutoffs) 
write.csv(ferro_grades_quantile, file =file.path(output_dir, "ferroptosis_grades_PAAD.csv"), row.names = TRUE)
Ferroptosis_grade<-ferro_grades_quantile

#Data processing
cox_data <- read.csv("my_data/cox_data.csv", header = TRUE, stringsAsFactors = FALSE)
cox_data_clean <- cox_data %>%
  mutate(
    histologic_grade = factor(histologic_grade, levels = c("G1", "G2", "G3", "G4", "GX")),
    residual_tumor = factor(residual_tumor, levels = c("R0", "R1", "R2", "RX")),
    stage_event = factor(stage_event, levels = c("Stage I", "Stage II", "Stage III", "Stage IV", "Stage X")),
    Ferroptosis_grade = factor(Ferroptosis_grade, levels=c("1","2","3","4","5","6","7","8","9","10")),
    status = ifelse(status == "Dead", 1, 0)
  ) %>%
  filter(
    histologic_grade != "GX",
    residual_tumor != "RX",
    stage_event != "Stage X"
  )%>%
  droplevels() 

#Only ferroptosis index as a standalone predictor
cox_data_clean$Ferroptosis_grade <- as.factor(cox_data_clean$Ferroptosis_grade)
dca_data <- decision_curve(
  status ~ Ferroptosis_grade, 
  data = cox_data_clean,
  family = binomial(link = "logit"),
  bootstraps = 500,
  confidence.intervals = 0.95)

pdf(file.path(output_dir, "DCA_Ferroptosisgrade.pdf"), width = 12, height = 8)
plot_decision_curve(dca_data)
dev.off()
plot_decision_curve(
  dca_data,
  col = c("#9370DB", "gray50", "gray70"),  
  legend.position = "none",
  confidence.intervals = TRUE,
  lty = c(1, 1, 1),  
  lwd = c(2, 2, 2)   
)

#Integrated prediction models, including tumor stage, ferroptosis index, treat all, treat none
dca_stage <- decision_curve(status ~ stage_event, data = cox_data_clean)
dca_ferro <- decision_curve(status ~ Ferroptosis_grade, data = cox_data_clean)

pdf(file.path(output_dir, "DCA_Integrated.pdf"), width = 12, height = 8)
par(mar = c(5, 6, 4, 2) + 0.1)  
plot(NULL, 
     xlim = c(0, 1), 
     ylim = c(-0.1, 1),  
     xlab = "High Risk Threshold", 
     ylab = "Standardized Net Benefit\n",
     cex.lab = 1.5,      
     cex.axis = 1.3,     
     cex.main = 1.8,    
     font.main = 2,      
     bty = "l",          
     las = 1)            
grid(col = "gray90", lty = 1, lwd = 0.5)
abline(h = 0, col = "#999999", lwd = 6, lty = 1) #None
lines(dca_stage$derived.data$thresholds, 
      dca_stage$derived.data$sNB.all,  
      col = "#666666", lwd = 6, lty = 1) #ALL
lines(dca_stage$derived.data$thresholds, 
      dca_stage$derived.data$sNB,  
      col = "#EACD76", lwd = 6, lty = 1) #Stage
lines(dca_ferro$derived.data$thresholds, 
      dca_ferro$derived.data$sNB,  
      col = "#745399", lwd = 6, lty = 1) #Ferroptosis Index
legend("topright",
       legend = c("Treat All", "Treat None", 
                  "Stage", "Ferroptosis"),
       col = c("#666666", "#999999", "#EACD76", "#745399"),
       lty = c(1, 1, 1, 1, 1), 
       lwd = 4,
       cex = 1.5,        
       bty = "n",        
       seg.len = 2,      
       x.intersp = 0.8,  
       y.intersp = 1.2)  
dev.off()











