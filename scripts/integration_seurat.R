#load library
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

#====== SPATIAL
#load all the files
list_value_files <- list.files(path=".", full.names=TRUE, pattern="_valuestp.txt", recursive = F)
names_value_files <- substr(list_value_files, 3,11)
#> names_value_files
#[1] "MIA01_PRE" "MIA08_PRE" "MIA09_PRE" "MIA10_PRE" "MIA11_EDT" "MIA11_PRE"
#[7] "MIA11_PRO"

#load all the file coordinates
list_xy_files <- list.files(path=".", full.names=TRUE, pattern="_XY.txt", recursive = F)
names_xy_files <- substr(list_xy_files, 3,14)
#> names_xy_files
#[1] "MIA01_PRE_XY" "MIA08_PRE_XY" "MIA09_PRE_XY" "MIA10_PRE_XY" "MIA11_EDT_XY"
#[6] "MIA11_PRE_XY" "MIA11_PRO_XY"

#create seurat object and downsample
#data_obj_downsample <- subset(data_obj, cells = sample(Cells(data_obj), 5000))
#assign(i, data_obj_downsample)
for(i in names_value_files){
data <- as.matrix(fread(paste0("/home/meltest/Data/CODEX/seurat_codex/23April2021/",i, "_valuestp.txt")), rownames=1)
data_obj <- CreateSeuratObject(counts = data, assay="Spatial")
assign(i, data_obj)
}

#Run per patient - mia11
#MIA11_EDT
MIA11_EDT@images$image =  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "MIA11EDT_",
    coordinates = MIA11_EDT_XY
  )


#MIA11_PRE
MIA11_PRE@images$image =  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "MIA11PRE_",
    coordinates = MIA11_PRE_XY
  )


#MIA11_PRO
MIA11_PRO@images$image =  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "MIA11PROG_",
    coordinates = MIA11_PRO_XY
  )

#merge only mia11 which has 3 time points
mia11.merge <- merge(MIA11_PRE, y = c(MIA11_PRO, MIA11_EDT), add.cell.ids = c("1MIA11PRE_CODEX", "2MIA11PRO_CODEX", "3MIA11EDT_CODEX"))

#mia11 is used to test out the integration.
#used the same normalisation as scRNAseq which is lognormalise method.
mia11.merge <- NormalizeData(mia11.merge, normalization.method = "LogNormalize", scale.factor = 10000, assay = "Spatial")
mia11.merge <- FindVariableFeatures(mia11.merge)
all.proteins <- rownames(mia11.merge)
mia11.merge <- ScaleData(mia11.merge, features = all.proteins)

#Run PCA
mia11.merge <- RunPCA(mia11.merge, approx=FALSE) #error due to the zero intensity of protein and read up on Seurat document as it is OK to ignore the error message

#Find clusters
mia11.merge <- FindNeighbors(mia11.merge, dims = 1:20) #can change to 1:30 for better res
mia11.merge <- FindClusters(mia11.merge, resolution = 0.4) 
mia11.merge <- RunUMAP(mia11.merge, dims = 1:20)

#Find markers
mia11.mergemarkers <- FindAllMarkers(mia11.merge, min.pct = 0.25, logfc.threshold = 0)


#=====scRNAseq
#load the single cell data but only MIA11
dat <- read.table("~/gene_counts/count_updated_MIA11.txt", row.names=1)
dat <- as.data.frame(dat)

#Read the metadata
mdata <- read.csv("~/metadata/melanoma_metadata_CQmia11.csv", row.names=1, sep=",")
mdata <- as.data.frame(mdata)

#Check that the names are the same in the data and metadata.
#Read count data and meta data not match because they are not in the same order. Run the codes below to match them up.
idx <- match(colnames(dat), rownames(mdata))
reordered_mdata <- mdata[idx,]
all(row.names(reordered_mdata) == colnames(dat))
 
#set up Seurat Object
mia11 <- CreateSeuratObject(counts = dat, min.cells =3, min.features = 200, meta.data = reordered_mdata, project = "singlecell_mia11", assay = "RNA")

#filter cells based on the above featurescatter graphs
mia11_filtered <- subset(mia11, subset = nGene < 7500 & percent.mt < 50)

#normalise data - default settings
mia11 <- NormalizeData(mia11_filtered, normalization.method = "LogNormalize", scale.factor = 10000)

#find variable features
mia11 <- FindVariableFeatures(mia11)

#Scale the data before draw tsne/umap
all.genes <- rownames(mia11)
mia11 <- ScaleData(mia11, features = all.genes)

#Run PCA
mia11 <- RunPCA(mia11, features = VariableFeatures(object = mia11))

#Cluster the cells
set.seed(1234)
mia11 <- FindNeighbors(mia11, dims = 1:20)
mia11 <- FindClusters(mia11, resolution = 0.4)
mia11 <- RunUMAP(mia11, dims = 1:20)

#Find all markers
mia11.markers <- FindAllMarkers(mia11, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#====== Anchoring
#Prior to plotting the UMAP, might need to downsample cells
mia11.merge_downsample <- subset(mia11.merge, cells = sample(Cells(mia11.merge), 15000)) 
mia11.merge_downsample_use <- as.vector(colnames(mia11.merge_downsample))


anchors <- FindTransferAnchors(reference = mia11, query = mia11.merge_downsample, normalization.method = "LogNormalize", npcs=20, k.anchor = 10, k.score = 50, approx=FALSE, k.filter=NA,  dims=1:20, max.features = NA)

predictions.assay <- TransferData(anchorset = anchors, refdata = mia11$seurat_clusters, k.weight = 50, n.trees = 50, weight.reduction = mia11.merge_downsample[["pca"]], dims = 1:20)

mia11.merge_downsample <- AddMetaData(mia11.merge_downsample, metadata = predictions.assay) #add metadata

predictions <- table(mia11.merge_downsample$seurat_clusters, mia11.merge_downsample$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

#Check and plot prediction score
correct <- length(which(mia11.merge$seurat_clusters == mia11.merge$predicted.id))
incorrect <- length(which(mia11.merge$seurat_clusters != mia11.merge$predicted.id))
data <- FetchData(mia11.merge_downsample, vars = c("prediction.score.max"))
pdf("predscore_codex_rna1.pdf")
ggplot(data, aes(prediction.score.max)) + geom_density(alpha = 0.5) + theme_cowplot() + xlab("Prediction Score")
dev.off()

#umap to show what clusters were anchored
png("umapanchor_codex_rna1.png", height=15, width=15, res=300, units="cm")
DimPlot(mia11.merge_downsample, group.by = "predicted.id", label = TRUE, raster=FALSE) + ggtitle("Predicted annotation after integration")
dev.off()

# Failed attempt for the below codes
#coembeded for visualisation
mia11$id <- 'singlecellRNAseq'
mia11.merge_downsample$id <- 'CODEX'
coembed <- merge(x = mia11, y = mia11.merge_downsample)

#coembed <- RunPCA(coembed)
#coembed <- RunUMAP(coembed, reduction = 'pca', dims = 1:20)
#DimPlot(coembed, group.by = "orig.ident", shuffle = TRUE)

=======
genes.use <- VariableFeatures(mia11)
protein.use <- VariableFeatures(mia11.merge) #this works due to small features
#ref_rna <- GetAssayData(mia11, assay="RNA", slot="data")[genes.use, ]

#imputation <- TransferData(anchorset = anchors, refdata = ref_rna, weight.reduction = mia11.merge[["pca"]], dims = 2:20)
#refquery
#add metadata
#mia11.merge[["RNA"]] <- imputation #add metadata
#mia11.merge <- AddMetaData(mia11.merge, metadata = imputation)

# run PCA and UMAP on this combined object, to visualize the co-embedding of both data sets
coembed <- ScaleData(coembed, features = protein.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = protein.use, verbose = FALSE, approx=FALSE)
coembed <- RunUMAP(coembed, dims = 1:20)

png("umapanchor_codexwithrna_DOWN2.png", height=15, width=30, res=300, units="cm")
DimPlot(coembed, group.by = c("id", "seurat_clusters"), raster=TRUE) #raster=TRUE 
dev.off()



