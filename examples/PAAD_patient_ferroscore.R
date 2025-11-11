######Calculate the ferroptosis score of PAAD patients
######################################################


library(FerroScore)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("edgeR", "stringr", "igraph", "ggraph", "STRINGdb",
                       "dplyr", "clusterProfiler", "org.Mm.eg.db", "org.Hs.eg.db",
                       "enrichplot", "ggplot2", "GOplot", "openxlsx"))

library(edgeR)
library(stringr)
library(igraph)
library(ggraph)
library(STRINGdb)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(GOplot)
library(openxlsx)



counts_file <- "my_data/PAAD_count_treated.csv"
meta_file <- "my_data/PAAD_meta_treated.csv"
output_dir <- "my_results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}



#Load pancreatic cancer data
counts_data <- read.csv(counts_file, row.names = 1)
meta_data <- read.csv(meta_file, stringsAsFactors = TRUE)


#Calculate ferroptosis scores for each patient
result_PAAD<- run_ferroptosis_analysis(
  counts_data = counts_data,
  meta_data = meta_data,
  control_group = "control",
  species = "human",
  ppi_file = NULL,
  output_dir = tempdir(),
  save_intermediate = FALSE,
  parallel = FALSE,
  n_cores = 2)


result_PAAD$summary
result_PAAD$merged_DEG
result_PAAD$results


#Save ferroptosis scores
ferroptosis_scores <- t(result_PAAD$summary)
write.csv(ferroptosis_scores, file = file.path(output_dir, "ferroptosis_scores_PAAD.csv"), row.names = TRUE)


#Save differential gene expression results
deg_list <- split(result_PAAD$merged_DEG, result_PAAD$merged_DEG$Sample)
names(deg_list)
wb <- createWorkbook()
for (sample_name in names(deg_list)) {
  addWorksheet(wb, sheetName = sample_name)
  writeData(wb, sheet = sample_name, x = deg_list[[sample_name]])
}
saveWorkbook(wb, file = file.path(output_dir, "diff_genes_PAAD.csv"), overwrite = TRUE)




