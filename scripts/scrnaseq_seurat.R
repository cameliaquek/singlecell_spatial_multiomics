#load library
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(SeuratWrappers)
library(patchwork)
library(Matrix)
library(celldex)
library(SingleR)

#Read the gene count matrix
dat <- read.table("~/gene_counts/count_updated_allpts1x.txt", row.names=1)
dat <- as.data.frame(dat)

#Read the metadata
mdata <- read.csv("~/melanoma_metadata_CQallpts.csv", row.names=1, sep=",")
mdata <- as.data.frame(mdata)

#Check that the names are the same in the data and metadata.
#Read count data and meta data not match because they are not in the same order. Run the codes below to match them up.
idx <- match(colnames(dat), rownames(mdata))
reordered_mdata <- mdata[idx,]
all(row.names(reordered_mdata) == colnames(dat))

#set up Seurat Object
immx <- CreateSeuratObject(counts = dat, min.cells =3, min.features = 200, meta.data = reordered_mdata, project = "singlecell_immx", assay = "RNA")

#Remove mia1b because this sample do not have matched samples from CODEX
immx.subset <- subset(immx, subset = meta.data@mia_ID != "MIA1b")
immx <- immx.subset

#If neeed to add CODEX ID into the meta data for other purposes - uncheck comments and run below example codes
#x <- immx@meta.data$codex_ID
#x <- gsub("mia11a","P1",x)
#immx@meta.data$codex_ID <- x

#Filter low quality cells
immx_filtered <- subset(immx, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 40)

#normalise data
immx <- NormalizeData(immx_filtered, normalization.method = "LogNormalize", scale.factor = 10000)

#find variable features
immx <- FindVariableFeatures(immx)

#Scale the data before draw tsne/umap
all.genes <- rownames(immx)
immx <- ScaleData(immx, features = all.genes)

#Run PCA and draw PCA heatmap/plot
immx <- RunPCA(immx, features = VariableFeatures(object = immx))

#Cluster the cells
set.seed(92021)
immx <- FindNeighbors(immx, dims = 1:20)
immx <- FindClusters(immx, resolution = 0.4)
immx <- RunUMAP(immx, dims = 1:20)

#Draw UMAP - uncheck comments if needed
#This parameter generates 18 clusters
#pdf("UMAPplots.pdf")
#DimPlot(immx, reduction = "umap")
#DimPlot(immx, reduction = "umap", label=TRUE)
#dev.off()

#Find all markers
immx.markers <- FindAllMarkers(immx, only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.15)
top10 <- immx.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#Draw heatmap to check on the clustering
DoHeatmap(immx, features = top10$gene) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))) + theme(text = element_text(size = 8))


#Run SingleR - After knowning that the Louvain clustering parameters fits well the study, we will begin characterising the cell subsets using singleR
#human refernece database
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se

immune.ref <- MonacoImmuneData()
immunestroma.ref <- BlueprintEncodeData()

#convert seurat object to singleR obj
objR <- as.SingleCellExperiment(immx)

#Run the databases
pred.hpca <- SingleR(test = objR, ref = hpca.se, assay.type.test=1, labels = hpca.se$label.main)
pred.immune.ref <- SingleR(test = objR, ref = immune.ref, assay.type.test=1, labels = immune.ref$label.fine)
pred.immunestroma.ref <- SingleR(test = objR, ref = immunestroma.ref, assay.type.test=1, labels =immunestroma.ref$label.fine)

#Check the results if needed
#head(sort(table(pred.hpca$labels), decreasing=TRUE), 10)
#head(sort(table(pred.immune.ref$labels), decreasing=TRUE), 10)
#head(sort(table(pred.immunestroma.ref$labels), decreasing=TRUE), 10)
#Alternatively, run the below example code to save to table
#write.table(pred.hpca, file="SingleR_results_hpca.txt", quote=FALSE, sep="\t", col.names=NA, row.names=TRUE)

