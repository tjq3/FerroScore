######Data processing and ferroptosis score calculation for AD patients
######################################################



library(FerroScore)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("edgeR", "stringr", "igraph", "ggraph", "STRINGdb",
                       "dplyr", "clusterProfiler",  "org.Hs.eg.db",
                       "enrichplot", "ggplot2", "GOplot", "openxlsx"))

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


output_dir <- "my_results"


#Load Alzheimer's disease data
#Data processing
counts_rawdata <- read.csv("my_data/AD_gene_all_counts_matrix.csv")
counts_rawdata_unique <- counts_rawdata %>%
  distinct(across(1), .keep_all = TRUE)  
write.csv(counts_rawdata_unique, file = file.path(output_dir, "AD_gene_all_counts_matrix_deduplication.csv"), row.names = TRUE)

#Load data
counts_data <- read.csv("my_data/AD_count_treated.csv", row.names = 1)
meta_data <- read.csv("my_data/AD_meta_treated.csv", stringsAsFactors = TRUE)



#Calculate ferroptosis scores for each patient
result_AD<- run_ferroptosis_analysis(
  counts_data = counts_data,
  meta_data = meta_data,
  control_group = "control",
  species = "human",
  ppi_file = NULL,
  output_dir = tempdir(),
  save_intermediate = FALSE,
  parallel = FALSE,
  n_cores = 2)


result_AD$summary
result_AD$merged_DEG
result_AD$results


#Save ferroptosis scores
ferroptosis_scores <- t(result_AD$summary)
write.csv(ferroptosis_scores, file = file.path(output_dir, "ferroptosis_scores_AD_treated.csv"), row.names = TRUE)


#Save differential gene expression results
deg_list <- split(result_AD$merged_DEG, result_AD$merged_DEG$Sample)
names(deg_list)
wb <- createWorkbook()
for (sample_name in names(deg_list)) {
  addWorksheet(wb, sheetName = sample_name)
  writeData(wb, sheet = sample_name, x = deg_list[[sample_name]])
}
saveWorkbook(wb, file = file.path(output_dir, "diff_genes_AD_treated.csv"), overwrite = TRUE)




