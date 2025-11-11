rm(list = ls())

getwd()
setwd("C:\\Users\\17865\\Documents\\FerroScore\\test")


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



counts_data <- read.csv("C:\\Users\\17865\\Documents\\FerroScore\\data\\GSE237928_dataAA.csv", row.names = 1)
meta_data <- read.csv("C:\\Users\\17865\\Documents\\FerroScore\\data\\GSE237928_metaAA.csv", stringsAsFactors = TRUE)

result_AA<- run_ferroptosis_analysis(
    counts_data = counts_data,
    meta_data = meta_data,
    control_group = "control",
    species = "mouse",
    ppi_file = NULL,
    output_dir = tempdir(),
    save_intermediate = FALSE,
    parallel = FALSE,
    n_cores = 2)

result_AA$summary
result_AA$merged_DEG
result_AA$results
