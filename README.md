

# FerroScore: A statistical approach for quantifying tumor-related ferroptosis based on omics data
FerroScore is an omics-based statistical method that quantifies ferroptosis by integrating activity scores from three core pathways—iron, glutathione, and lipid metabolism—into a continuous score. This framework enables cross-resolution assessment of ferroptosis and provides mechanistic insights into tumors, immune processes, and neurodegenerative diseases, thus having potential applications in targeted therapy and drug discovery.

# Overview
FerroScore is an algorithm that quantifies the ferroptosis activity based on transcriptomic data (RNA-seq or scRNA-seq). It enables quantitative assessment of ferroptosis in individual samples or single cells by systematically integrating key molecular features and network topology related to ferroptosis.
![image](https://github.com/tjq3/FerroScore/blob/main/overview_FerroScore.jpg)
FerroScore converts transcriptomic data into quantitative ferroptosis metrics. The pipeline takes RNA-seq or scRNA-seq data and associated labels as input (“**Input**”). It constructs a PPI network from DEGs, refines it into a core ferroptosis network via module identification, functional annotation, and merging, and uses topological features to generate a continuous activity score and a stratified index (range: 1-10) (“**FerroScore**”).  These outputs enable diverse downstream analyses, including association with clinical survival, identification of high-activity cell subpopulations, and inference of autocrine and paracrine regulatory networks (“**Output**”).

# Code Structure
The respository mainly consists of the following modules:

- `inst` provides PPI files for both mouse (Mus musculus) and human (Homo sapiens) species.
- `calculate_ferroptosis.R` is used to calculate the ferroptosis score for individual samples (tissues or cells). This function performs differential expression analysis, constructs the corresponding PPI network, conducts module identification and functional annotation, extracts the ferroptosis-related network from it, calculates the topological features of this network, and thereby derives the ferroptosis score.
- `run_ferroptosis_analysis.R` is used to run ferroptosis analysis for multiple samples (tissues or cells).
- `data` contains the dataset used for benchmarking FerroScore on the experimental dataset, using AA as an example.
- `test` contains the code used for benchmarking FerroScore on the experimental dataset, using AA as an example.
- `examples` illustrates the application of FerroScore across several domains: specifically in pancreatic cancer at the patient, cellular, and immune levels, as well as in Alzheimer's disease. The provided content encompasses the data preprocessing pipelines for diverse datasets and the ensuing downstream analyses.

# Installation
## Install from Github using devtools
```
devtools :: install_github("tjq3/FerroScore")
```
## Install from R source codes
Download source codes and type in R
```
install.packages(path_to_file, type = 'source', rep = NULL) # The path_to_file would represent the full path and file name
```
## Dependencies
FerroScore relies on several libraries. Please install any other dependencies if they are not installed automatically.

- edgeR
- stringr
- igraph
- ggraph
- STRINGdb
- dplyr
- clusterProfiler
- org.Mm.eg.db
- org.Hs.eg.db
- enrichplot
- ggplot2
- GOplot
- openxlsx

# Datasets
All datasets can be downloaded from their respective sources:

- The mouse MEF cell line RNA-seq data used in this study are available from the Gene Expression Omnibus (GEO; [https://www.ncbi.nlm.nih.gov/geo/](https://www.ncbi.nlm.nih.gov/geo/)) under accession number GSE237928.
- The human A549 cell line RNA-seq data can be obtained from the Sequence Read Archive (SRA; [https://www.ncbi.nlm.nih.gov/sra](https://www.ncbi.nlm.nih.gov/sra)) under accession number PRJNA979805.
- Bulk RNA-seq data and corresponding clinical information for pancreatic cancer are publicly accessible via The Cancer Genome Atlas (TCGA; [https://portal.gdc.cancer.gov/](https://portal.gdc.cancer.gov/)).
- scRNA-seq data of immune and non-immune cells from pancreatic tumor tissues are downloaded from the GEO database under accession numbers GSE235449 and GSE194247, respectively.
- RNA-seq and clinical data for Alzheimer’s disease are sourced from The Mayo RNAseq Study (MayoRNAseq) project within the AD Knowledge Portal ([https://adknowledgeportal.org](https://adknowledgeportal.org/)).
