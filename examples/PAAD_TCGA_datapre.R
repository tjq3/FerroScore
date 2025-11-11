#####TCGA-PAAD RNA-seq data preprocess
#############################################


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rjson")
BiocManager::install("tidyverse")

library(rjson)
library(tidyverse)



################################ Count
#Read metadata file
json<-jsonlite::fromJSON("metadata.cart.2025-03-11.json")
view(json)

#Get sample names and file names
sample_id<-sapply(json$associated_entities,function(x){x[,1]})
sample_id[1:10]
file_sample<-data.frame(sample_id,file_name=json$file_name)
view(file_sample)

#Get count file locations
count_file<-list.files('gdc_download_20250311_031921.130048/',
                       pattern = '*.tsv',recursive = TRUE)
count_file[1:10]

#Get each file name
count_file_name<-strsplit(count_file,split = '/')
count_file_name<-sapply(count_file_name, function(x){x[2]})
count_file_name[1:10]

#Create empty data frame
matrix=data.frame(matrix(nrow = 60660,ncol = 0))

#Read and merge files one by one
for (i in 1:length(count_file)) {
  path=paste0('gdc_download_20250311_031921.130048//',count_file[i])
  data<-read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[2,]
  data<-data[-c(1:6),]
  data<-data[3]  ## 3 represents unnormalized count; 4 represents once normalized; 5 represents twice normalized; 6 represents TPM; 7 represents FPKM; 8 represents FPKM_uq
  colnames(data)<-file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix<-cbind(matrix,data)
}

#Convert to gene symbols
path=paste0('gdc_download_20250311_031921.130048//',count_file[1])
data<-as.matrix(read.delim(path,fill = TRUE,header = FALSE,row.names = 1))
gene_name<-data[-c(1:6),1]
gene_name[1:10]
matrix0<-cbind(gene_name,matrix)
gene_type<-data[-c(1:6),2]
gene_type[1:10]
matrix0<-cbind(gene_type,matrix0)

#Remove duplicate genes, keep the result with maximum expression
matrix0<-aggregate(.~gene_name,data=matrix0,max)
table(gene_name)

#Keep only mRNA
matrix0<-subset(x=matrix0,gene_type=="protein_coding")
table(gene_type)

#Set gene_name column as row names and convert to export format
rownames(matrix0)<-matrix0[,1]
matrix0<-matrix0[,-c(1,2)]
matrix1=data.frame(ID=rownames(matrix0),matrix0)
colnames(matrix1)=gsub('[.]','-',colnames(matrix1))

write.csv(matrix1,file="TCGA_PAAD_count.csv",quote = F)








################################## TPM
json<-jsonlite::fromJSON("metadata.cart.2025-03-11.json")
view(json)

sample_id<-sapply(json$associated_entities,function(x){x[,1]})
sample_id[1:10]
file_sample<-data.frame(sample_id,file_name=json$file_name)
view(file_sample)

count_file<-list.files('gdc_download_20250311_031921.130048/',
                       pattern = '*.tsv',recursive = TRUE)
count_file[1:10]

count_file_name<-strsplit(count_file,split = '/')
count_file_name<-sapply(count_file_name, function(x){x[2]})
count_file_name[1:10]

matrix=data.frame(matrix(nrow = 60660,ncol = 0))

for (i in 1:length(count_file)) {
  path=paste0('gdc_download_20250311_031921.130048//',count_file[i])
  data<-read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[2,]
  data<-data[-c(1:6),]
  data<-data[6]  ## 3 represents unnormalized count; 4 represents once normalized; 5 represents twice normalized; 6 represents TPM; 7 represents FPKM; 8 represents FPKM_uq
  colnames(data)<-file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix<-cbind(matrix,data)
}

path=paste0('gdc_download_20250311_031921.130048//',count_file[1])
data<-as.matrix(read.delim(path,fill = TRUE,header = FALSE,row.names = 1))
gene_name<-data[-c(1:6),1]
gene_name[1:10]
matrix0<-cbind(gene_name,matrix)
gene_type<-data[-c(1:6),2]
gene_type[1:10]
matrix0<-cbind(gene_type,matrix0)

matrix0<-aggregate(.~gene_name,data=matrix0,max)
table(gene_name)

matrix0<-subset(x=matrix0,gene_type=="protein_coding")
table(gene_type)

rownames(matrix0)<-matrix0[,1]
matrix0<-matrix0[,-c(1,2)]
matrix1=data.frame(ID=rownames(matrix0),matrix0)
colnames(matrix1)=gsub('[.]','-',colnames(matrix1))

write.csv(matrix1,file="TCGA_PAAD_TPM.csv",quote = F)




















