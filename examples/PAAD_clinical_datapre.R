######clinical-PAAD RNA-seq data preprocess
#############################################


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("XML")
BiocManager::install("methods")

library(XML)
library(methods)



################################ clinical
dir <- "gdc_download_20250311_034900.311087"

#Get the location of each sample
all_files=list.files(path = dir,pattern = '*.xml$',recursive = T)
all_files[1:10]

cl=lapply(all_files,function(x){
  #Read XML file
  result<-xmlParse(file = file.path(dir,x))
  #Get root element
  rootnode<-xmlRoot(result)
  #Convert to data frame
  xmldataframe<-xmlToDataFrame(rootnode[2])
  return(t(xmldataframe))
})

#Merge each list
clinical<-t(do.call(cbind,cl))

write.csv(clinical,file="TCGA_PAAD_clinical.csv",quote = F)


