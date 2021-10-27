#load library
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)

#Akoya team has CODEX data and thus below is an example file and R examples used to run the analysis
mia11pre_data <- as.matrix(fread("Run_001_MIA11_Precell_matrix1.txt"),rownames=1)
mia11pre_coord <- read.table("Run_001_MIA11_PREcell_xy.txt", header=TRUE)
mia11pre_coord.df <- data.frame(x=mia11pre_coord$x, y=mia11pre_coord$y, stringsAsFactors=FALSE)
rownames(mia11pre_coord.df) <- mia11pre_coord$cell_id

#In order for Seurat to read spatial data, below codes are needed for Seurat to read spatial information
mia11pre <- CreateSeuratObject(counts = mia11pre_data, assay="Spatial")
mia11pre@images$image =  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "MIA11PRE_",
    coordinates = mia11pre_coord.df
  )

#downsample cells to 3000
#Choices of normalisation including SCTransform and log transformation. SCTransform used here.
mia11pre <- SCTransform(mia11pre, assay = "Spatial", ncells=3000)

#Find clusters
mia11pre <- RunPCA(mia11pre, approx=FALSE) #error due to the zero intensity of protein and read up on Seurat document as it is OK to ignore the error message
mia11pre <- FindNeighbors(mia11pre, dims = 1:20)
mia11pre <- FindClusters(mia11pre, resolution = 0.4) 
mia11pre <- RunUMAP(mia11pre, dims = 1:20)
