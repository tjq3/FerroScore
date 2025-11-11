######Clinical Characteristics and Mechanisms of Ferroptosis in AD Patients
######################################################



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("tidyverse", "readr","ggplot2","mgcv","rstatix",
                       "PMCMRplus","vcd","ggpubr","coin","reshape2",
                       "pheatmap","clusterProfiler","org.Hs.eg.db",
                       "ggrepel","tidytext","dplyr","patchwork","tibble"))

library(tidyverse)
library(readr)
library(ggplot2)
library(mgcv)
library(rstatix)
library(PMCMRplus)
library(vcd)       
library(ggpubr)
library(coin)
library(reshape2)
library(pheatmap)     
library(clusterProfiler)  
library(org.Hs.eg.db) 
library(ggrepel)      
library(tidytext)  
library(dplyr)
library(patchwork)
library(tibble)


output_dir <- "my_results/AD"


####Ferroptosis Score & Death Age
ferro_data <- read_csv("my_data/ferroptosis_scores_AD_treated.csv")
ferro_scores <- ferro_data[4,]
ferro_scores_transposed <- ferro_scores %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "patient_id")
ferro_scores<-ferro_scores_transposed[-1,]
colnames(ferro_scores)[2] <- "Ferroptosis_score"

clinical_data <- read_csv("my_data/AD_clinical.csv")
death_age <- clinical_data %>% 
  select(patient_id = 1, death_age = 2)

merged_data_deathage <- inner_join(death_age, ferro_scores, by = "patient_id")
write.csv(merged_data_deathage, file = file.path(output_dir, "merged_data_AD_deathage.csv"), row.names = TRUE)

pdf(file.path(output_dir, "FerroScores&DeathAge_empiricalcurve.pdf"), width = 9,height = 6)
ggplot(merged_data_deathage, aes(x = as.numeric(Ferroptosis_score), y = death_age)) +
  geom_point(alpha = 1, size = 3, color = "#154996") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), 
              color = "#7291c0", se = TRUE, fill = "#b8c8df",
              linetype = "solid", level = 0.95) +
  labs(x = "Ferroptosis Score", y = "Death Age") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        panel.border = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.8),  
        axis.text = element_text(color = "black", size = 13),
        axis.title = element_text(face = "bold", size = 18),
        plot.margin = margin(15, 15, 15, 15)
  )
dev.off()



#####Ferroptosis Score & APOE
apoe <- clinical_data %>% 
  select(patient_id = 1, apoe = 4)

merged_data_apoe <- inner_join(apoe, ferro_scores, by = "patient_id")
write.csv(merged_data_apoe, file = file.path(output_dir, "merged_data_AD_apoe.csv"), row.names = TRUE)

analysis_data_as <- merged_data_apoe %>%
  mutate(
    apoe = factor(apoe, 
                  levels = c("23", "33", "34", "44"),
                  labels = c("E2E3", "E3E3", "E3E4", "E4E4"),
                  ordered = TRUE),
    Ferroptosis_score = as.numeric(Ferroptosis_score)
  ) %>%
  filter(!is.na(apoe), !is.na(Ferroptosis_score))

pdf(file.path(output_dir, "FerroScores&APOE_box.pdf"),width = 9,height = 6)
grade_colors <- c("E2E3" = "#FFCCCC",  
                  "E3E3" = "#FFE6CC",   
                  "E3E4" = "#FFFFCC",    
                  "E4E4" = "#CCFFCC")    
ggplot(analysis_data_as, aes(x = apoe, y = Ferroptosis_score, fill = apoe)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0.8 ) +
  geom_jitter(width = 0.15, alpha = 0.55, size = 3, shape = 21, color = "black",
              stroke = 0.55, aes(color = apoe), show.legend = FALSE) +
  scale_fill_manual( values = grade_colors) +
  scale_color_manual(values = grade_colors)+
  labs(x = "APOE Genotype", 
       y = "Ferroptosis Score") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        panel.border = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.8),  
        axis.text = element_text(color = "black", size = 13),
        axis.title = element_text(face = "bold", size = 18),
        plot.margin = margin(15, 15, 15, 15),
        legend.position = "none"
  )
dev.off()



#####Ferroptosis Score & Braak
braak <- clinical_data %>% 
  select(patient_id = 1, braak = 5)

merged_data_braak <- inner_join(braak, ferro_scores, by = "patient_id")
write.csv(merged_data_braak, file = file.path(output_dir, "merged_data_AD_braak.csv"), row.names = TRUE)

analysis_data_bs <- merged_data_braak %>%
  mutate(
    braak = factor(braak, 
                   levels = c("4.5", "5", "5.5", "6"),
                   labels = c("Braak IV.5", "Braak V", "Braak V.5", "Braak VI"),
                   ordered = TRUE),
    Ferroptosis_score = as.numeric(Ferroptosis_score)
  ) %>%
  filter(!is.na(braak), !is.na(Ferroptosis_score))

pdf(file.path(output_dir, "FerroScores&Braak_box.pdf"),width = 9,height = 6)
grade_colors <- c("Braak IV.5" = "#FFCCCC",  
                  "Braak V" = "#FFE6CC",   
                  "Braak V.5" = "#FFFFCC",    
                  "Braak VI" = "#CCFFCC")    
ggplot(analysis_data_bs, aes(x = braak, y = Ferroptosis_score, fill = braak)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0.8 ) +
  geom_jitter(width = 0.15, alpha = 0.55, size = 3, shape = 21, color = "black",
              stroke = 0.55, aes(color = braak), show.legend = FALSE) +
  scale_fill_manual( values = grade_colors) +
  scale_color_manual(values = grade_colors)+
  labs(x = "Braak Stage", 
       y = "Ferroptosis Score") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        panel.border = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.8),  
        axis.text = element_text(color = "black", size = 13),
        axis.title = element_text(face = "bold", size = 18),
        plot.margin = margin(15, 15, 15, 15),
        legend.position = "none"
  )
dev.off()



#####Ferroptosis Score & Thal
thal <- clinical_data %>% 
  select(patient_id = 1, thal = 6)

merged_data_thal <- inner_join(thal, ferro_scores, by = "patient_id")
write.csv(merged_data_thal, file = file.path(output_dir, "merged_data_AD_thal.csv"), row.names = TRUE)

analysis_data_ts <- merged_data_thal %>%
  mutate(
    thal = factor(thal, 
                  levels = c("3", "4", "5"),
                  labels = c("Thal III","Thal IV","Thal V"),
                  ordered = TRUE),
    Ferroptosis_score = as.numeric(Ferroptosis_score)
  ) %>%
  filter(!is.na(thal), !is.na(Ferroptosis_score))

pdf(file.path(output_dir, "FerroScores&Thal_box.pdf"),width = 9,height = 6)
grade_colors <- c("Thal III" = "#FFCCCC",  
                  "Thal IV" = "#FFE6CC",   
                  "Thal V" = "#FFFFCC")    
ggplot(analysis_data_ts, aes(x = thal, y = Ferroptosis_score, fill = thal)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0.8 ) +
  geom_jitter(width = 0.15, alpha = 0.55, size = 3, shape = 21, color = "black",
              stroke = 0.55, aes(color = thal), show.legend = FALSE) +
  scale_fill_manual( values = grade_colors) +
  scale_color_manual(values = grade_colors)+
  labs(x = "Thal Phase", 
       y = "Ferroptosis Score") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        panel.border = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.8),  
        axis.text = element_text(color = "black", size = 13),
        axis.title = element_text(face = "bold", size = 18),
        plot.margin = margin(15, 15, 15, 15),
        legend.position = "none"
  )
dev.off()



#Ferroptosis characteristics of patients across different Braak
#DEGs in patients across different Braak
braak4.5_samples <- read_csv("my_data/braak4.5_id.csv")  
braak4.5_deg_files <- list.files("my_data/braak4.5_DEG", pattern = "\\.csv$", full.names = TRUE)
braak4.5_deg_all <- map_dfr(braak4.5_deg_files, ~ {
  read_csv(.x) %>%
    mutate(patient_id = str_extract(basename(.x), "X\\d+_(CER|TCX)")) 
})
braak4.5_deg_sig <- braak4.5_deg_all %>%
  filter(change != "Not Change")
write.csv(braak4.5_deg_sig, file = file.path(output_dir, "deg_sig_braak4.5.csv"), quote = F)


braak5_samples <- read_csv("my_data/braak5_id.csv")  
braak5_deg_files <- list.files("my_data/braak5_DEG", pattern = "\\.csv$", full.names = TRUE)
braak5_deg_all <- map_dfr(braak5_deg_files, ~ {
  read_csv(.x) %>%
    mutate(patient_id = str_extract(basename(.x), "X\\d+_(CER|TCX)")) 
})
braak5_deg_sig <- braak5_deg_all %>%
  filter(change != "Not Change")
write.csv(braak5_deg_sig, file = file.path(output_dir, "deg_sig_braak5.csv"), quote = F)


braak5.5_samples <- read_csv("my_data/braak5.5_id.csv")  
braak5.5_deg_files <- list.files("my_data/braak5.5_DEG", pattern = "\\.csv$", full.names = TRUE)
braak5.5_deg_all <- map_dfr(braak5.5_deg_files, ~ {
  read_csv(.x) %>%
    mutate(patient_id = str_extract(basename(.x), "X\\d+_(CER|TCX)")) 
})
braak5.5_deg_sig <- braak5.5_deg_all %>%
  filter(change != "Not Change")
write.csv(braak5.5_deg_sig, file = file.path(output_dir, "deg_sig_braak5.5.csv"), quote = F)


braak6_samples <- read_csv("my_data/braak6_id.csv")  
braak6_deg_files <- list.files("my_data/braak6_DEG", pattern = "\\.csv$", full.names = TRUE)
braak6_deg_all <- map_dfr(braak6_deg_files, ~ {
  read_csv(.x) %>%
    mutate(patient_id = str_extract(basename(.x), "X\\d+_(CER|TCX)")) 
})
braak6_deg_sig <- braak6_deg_all %>%
  filter(change != "Not Change")
write.csv(braak6_deg_sig, file = file.path(output_dir, "deg_sig_braak6.csv"), quote = F)



#Summarize all DEG information from all groups
deg_summary <- bind_rows(
  braak4.5_deg_sig %>% mutate(group = "Braak IV.5"),
  braak5_deg_sig %>% mutate(group = "Braak V"),
  braak5.5_deg_sig %>% mutate(group = "Braak V.5"),
  braak6_deg_sig %>% mutate(group = "Braak VI")
) %>%
  count(group, change, name = "count") %>%
  complete(group, change, fill = list(count = 0)) 
write.csv(deg_summary, file = file.path(output_dir, "deg_summary_braak.csv"), quote = F)

#Extract consistently high-frequency DEGs (consistent change in at least 50% of patients)
get_consensus_genes <- function(deg_data, min_patients_ratio = 0.5) {
  deg_data %>%
    group_by(Gene, change) %>%
    summarise(n_patients = n_distinct(patient_id), .groups = "drop") %>%
    filter(n_patients >= max(n_patients) * min_patients_ratio) %>%
    arrange(desc(n_patients))
}
braak4.5_consensus <- get_consensus_genes(braak4.5_deg_sig)
braak5_consensus <- get_consensus_genes(braak5_deg_sig)
braak5.5_consensus <- get_consensus_genes(braak5.5_deg_sig)
braak6_consensus <- get_consensus_genes(braak6_deg_sig)

#Combine high-frequency consensus genes from all groups
all_consensus <- bind_rows(
  braak4.5_consensus %>% mutate(group = "Braak IV.5"),
  braak5_consensus %>% mutate(group = "Braak V"),
  braak5.5_consensus %>% mutate(group = "Braak V.5"),
  braak6_consensus %>% mutate(group = "Braak VI")
) 
write.csv(all_consensus, file = file.path(output_dir, "consensusgene_all_braak.csv"), quote = F)

#Extract top 20 high-frequency consensus genes from each group
top_genes <- bind_rows(
  braak4.5_consensus %>% head(20) %>% mutate(group = "Braak IV.5"),
  braak5_consensus %>% head(20) %>% mutate(group = "Braak V"),
  braak5.5_consensus %>% head(20) %>% mutate(group = "Braak V.5"),
  braak6_consensus %>% head(20) %>% mutate(group = "Braak VI")
) %>% distinct(Gene, change, n_patients, group)
write.csv(top_genes, file = file.path(output_dir, "consensusgene_top_braak.csv"), quote = F)

#Create grouped bar chart
pdf(file.path(output_dir, "top_consensusgene_barchat_braak.pdf"), width = 14, height = 10)
ggplot(top_genes, aes(x = reorder_within(Gene, n_patients, group), 
                      y = n_patients, fill = change)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +  
  scale_fill_manual(values = c("Up" = "#e12216", "Down" = "#154996"),  
                    name = "Expression Change",  
                    labels = c("Down", "Up")) +  
  scale_x_reordered() +
  facet_wrap(~ group, scales = "free_y", ncol = 2) +  
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



###Heatmap of top DEGs
#Load raw expression matrix
expr_matrix <- read_csv("my_data/AD_gene_all_counts_matrix.csv")
valid_matrix<-read_csv("my_data/ferroptosis_scores_AD_treated.csv")

#Filter raw expression matrix for high-frequency consensus genes
expr_matrix <- as_tibble(expr_matrix)
heatmap_data <- expr_matrix %>%
  dplyr::select(Gene, any_of(colnames(valid_matrix))) %>%
  filter(Gene %in% top_genes$Gene) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

#Add annotation information (sample-group mapping)
annotation_col <- data.frame(
  Group = ifelse(colnames(heatmap_data) %in% braak4.5_samples$patient_id, "Braak IV.5",
                 ifelse(colnames(heatmap_data) %in% braak5_samples$patient_id, "Braak V",
                        ifelse(colnames(heatmap_data) %in% braak5.5_samples$patient_id, "Braak V.5", "Braak VI"))))
rownames(annotation_col) <- colnames(heatmap_data)

#Create heatmap
pdf(file.path(output_dir, "top_consensusgene_heatmap_braak.pdf"), width = 14, height = 10)
heatmap_colors <- colorRampPalette(c("#154996", "#f7f7f7", "#e12216"))(100)
annotation_colors <- list(
  Group = c("Braak IV.5" = "#FFCCCC", "Braak V" = "#FFE6CC", "Braak V.5" = "#FFFFCC", "Braak VI" = "#CCFFCC"))
pheatmap(
  heatmap_data,
  scale = "row",
  color = heatmap_colors,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors, 
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



###GO Enrichment Analysis
braak4.5_go_genes <- braak4.5_deg_sig %>% pull(Gene)
braak4.5_go_ids <- bitr(braak4.5_go_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
braak4.5_go_enrich <- enrichGO(
  gene = braak4.5_go_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
braak4.5_transGO <- setReadable(braak4.5_go_enrich, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(braak4.5_transGO, file = file.path(output_dir, "GO_braak4.5.csv"), quote = F)

braak5_go_genes <- braak5_deg_sig %>% pull(Gene)
braak5_go_ids <- bitr(braak5_go_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
braak5_go_enrich <- enrichGO(
  gene = braak5_go_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
braak5_transGO <- setReadable(braak5_go_enrich, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(braak5_transGO, file = file.path(output_dir, "GO_braak5.csv"), quote = F)

braak5.5_go_genes <- braak5.5_deg_sig %>% pull(Gene)
braak5.5_go_ids <- bitr(braak5.5_go_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
braak5.5_go_enrich <- enrichGO(
  gene = braak5.5_go_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
braak5.5_transGO <- setReadable(braak5.5_go_enrich, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(braak5.5_transGO, file = file.path(output_dir, "GO_braak5.5.csv"), quote = F)

braak6_go_genes <- braak6_deg_sig %>% pull(Gene)
braak6_go_ids <- bitr(braak6_go_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
braak6_go_enrich <- enrichGO(
  gene = braak6_go_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
braak6_transGO <- setReadable(braak6_go_enrich, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(braak6_transGO, file = file.path(output_dir, "GO_braak6.csv"), quote = F)

#Draw GO plots
p1 <- barplot(braak4.5_transGO, drop = TRUE, showCategory = 15,label_format = 100) + 
  ggtitle("Braak IV.5") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 21),  
        axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        legend.text = element_text(size = 16),  
        legend.title = element_text(size = 18), 
        panel.spacing = unit(1, "lines"),  
        plot.margin = margin(15, 15, 15, 15))
p2 <- barplot(braak5_transGO, drop = TRUE, showCategory = 15,label_format = 100) + 
  ggtitle("Braak V") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 21),  
        axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        legend.text = element_text(size = 16),  
        legend.title = element_text(size = 18), 
        panel.spacing = unit(1, "lines"),  
        plot.margin = margin(15, 15, 15, 15))
p3 <- barplot(braak5.5_transGO, drop = TRUE, showCategory = 15,label_format = 100) + 
  ggtitle("Braak V.5") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 21),  
        axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        legend.text = element_text(size = 16),  
        legend.title = element_text(size = 18), 
        panel.spacing = unit(1, "lines"),  
        plot.margin = margin(15, 15, 15, 15))
p4 <- barplot(braak6_transGO, drop = TRUE, showCategory = 15,label_format = 100) + 
  ggtitle("Braak VI") +
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

pdf(file.path(output_dir, "GO_braak.pdf"), width = 24, height = 16)
print(combined_plot)
dev.off()



###KEGG Enrichment Analysis
braak4.5_kegg <- enrichKEGG(gene = braak4.5_go_ids$ENTREZID,
                            organism = "hsa",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
braak4.5_kk <- setReadable(braak4.5_kegg, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(braak4.5_kk, file = file.path(output_dir, "KEGG_braak4.5.csv"), quote = F)

braak5_kegg <- enrichKEGG(gene = braak5_go_ids$ENTREZID,
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
braak5_kk <- setReadable(braak5_kegg, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(braak5_kk, file = file.path(output_dir, "KEGG_braak5.csv"), quote = F)

braak5.5_kegg <- enrichKEGG(gene = braak5.5_go_ids$ENTREZID,
                            organism = "hsa",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
braak5.5_kk <- setReadable(braak5.5_kegg, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(braak5.5_kk, file = file.path(output_dir, "KEGG_braak5.5.csv"), quote = F)

braak6_kegg <- enrichKEGG(gene = braak6_go_ids$ENTREZID,
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)
braak6_kk <- setReadable(braak6_kegg, OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')
write.csv(braak6_kk, file = file.path(output_dir, "KEGG_braak6.csv"), quote = F)

#Draw KEGG plots
p1 <- dotplot(braak4.5_kegg, showCategory = 15,label_format = 100) + 
  ggtitle("Braak IV.5") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 21),  
        axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        legend.text = element_text(size = 16),  
        legend.title = element_text(size = 18), 
        panel.spacing = unit(1, "lines"),  
        plot.margin = margin(15, 15, 15, 15))
p2 <- dotplot(braak5_kegg, showCategory = 15,label_format = 100) + 
  ggtitle("Braak V") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 21),  
        axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        legend.text = element_text(size = 16),  
        legend.title = element_text(size = 18), 
        panel.spacing = unit(1, "lines"),  
        plot.margin = margin(15, 15, 15, 15))
p3 <- dotplot(braak5.5_kegg, showCategory = 15,label_format = 100) + 
  ggtitle("Braak V.5") +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 21),  
        axis.text.y = element_text(color = "black", size = 18),
        axis.text.x = element_text(color = "black", size = 18),
        legend.text = element_text(size = 16),  
        legend.title = element_text(size = 18), 
        panel.spacing = unit(1, "lines"),  
        plot.margin = margin(15, 15, 15, 15))
p4 <- dotplot(braak6_kegg, showCategory = 15,label_format = 100) + 
  ggtitle("Braak VI") +
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

pdf(file.path(output_dir, "KEGG_braak.pdf"), width = 24, height = 16)
print(combined_plot)
dev.off()







