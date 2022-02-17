#Load library
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

#Analysis for all patient
all <- readRDS("/Users/camelia/Projects/immx_evolution/CITEseq/MIA_ALL_obj.rds")

#Filter negative and doublet from analysis
all_filter <- subset(all, subset = sample_ID != "Doublet" | sample_ID != "Negative")

#Sample ID should contains the names of all samples if mouse/genative/doublet are removed.
#> unique(all_filter@meta.data$sample_ID)
#[1] "MELCAP01_MIA01_LN_Pre"        "MELCAP01_MIA02_Subcut_TBC"   
#[3] "MELCAP02_MIA11_Subcut_Pre-1"  "MELCAP02_MIA01_LN_Prog"      
#[5] "MELCAP03_MIA11-Subcut_EDT-3"  "MELCAP03_MIA09_LN_Pre"       
#[7] "MELCAP04_MIA10_Subcut_Pre"    "MELCAP04_MIA11-Subcut_Prog-2"

#/////////////////redo anchoring on IMMUNE cells
# Remove mouse/rat labels
all_cite <- GetAssayData(object = all, assay = "ADT", slot = "data")
all_rna <- GetAssayData(object = all, assay = "RNA", slot = "data")
all.rna1 <- all.rna[!grepl("mm10", row.names(all.rna)),]
all.adt1 <- all.adt[!grepl("CITE-Mouse|CITE-Rat|CITE-Arm.ham", row.names(all.adt)),]
all.rna <- all.rna1
all.adt <- all.adt1

# Check names - since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(all.rna), colnames(all.adt))

#metadata
all_meta <- as.data.frame(all@meta.data)
idx <- match(colnames(all.rna), rownames(all_meta))
reordered_mdata <- all_meta[idx,]
all(row.names(reordered_mdata) == colnames(all.rna))

# creates a Seurat object based on the scRNA-seq data
ptall <- CreateSeuratObject(counts = all.rna, meta.data = reordered_mdata)
adt_assay <- CreateAssayObject(counts = all.adt)
ptall[["ADT"]] <- adt_assay

# Validate that the object now contains multiple assays
Assays(ptall)
#> Assays(ptall)
#[1] "RNA" "ADT"

#Running WNN
DefaultAssay(ptall) <- 'RNA'
ptall <- NormalizeData(ptall) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(ptall) <- 'ADT'
VariableFeatures(ptall) <- rownames(ptall[["ADT"]])

ptall <- NormalizeData(ptall, normalization.method = "CLR", margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = "apca")

ptall <- FindMultiModalNeighbors(
  ptall, reduction.list = list("pca", "apca"), 
  dims.list = list(1:25, 1:18), modality.weight.name = "RNA.weight")
  
ptall <- RunUMAP(ptall, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
ptall <- FindClusters(ptall, graph.name = "wsnn", algorithm = 1, resolution = 2)
ptall <- RunUMAP(ptall, reduction = 'pca', dims = 1:30, assay = 'RNA', 
              reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
ptall <- RunUMAP(ptall, reduction = 'apca', dims = 1:18, assay = 'ADT', 
              reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

pdf("UMAP-WNN_all_allcells.pdf")
DimPlot(ptall, reduction = "wnn.umap", label = FALSE)
DimPlot(ptall, reduction = "wnn.umap", label = TRUE)
DimPlot(ptall, reduction = "wnn.umap", group.by = "sampleID", label = FALSE)
DimPlot(ptall, reduction = "wnn.umap", group.by = "sampleID", label = TRUE)
DimPlot(ptall, reduction = "wnn.umap", split.by = "sampleID", label = FALSE)
DimPlot(ptall, reduction = "wnn.umap", split.by = "sampleID", label = TRUE)
dev.off()
p3 <- DimPlot(ptall, reduction = 'rna.umap', label = TRUE, 
              repel = TRUE, label.size = 2.5) + ggtitle("RNA_scRNAseq")
p4 <- DimPlot(ptall, reduction = 'adt.umap', label = TRUE, 
              repel = TRUE, label.size = 2.5) + ggtitle("Protein_CITEseq")
p5 <- DimPlot(ptall, reduction = "wnn.umap",
              label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
png("UMAP-all_allcells25-15-0.8.png", width=35, height=15, units="cm", res=500)
p3 + p4 + p5
dev.off()


DefaultAssay(ptall) <- "ADT"
p6 <- FeaturePlot(ptall, features = c("CITE-CD3", "CITE-CD8A", "CITE-PD1", "CITE-CTLA4"),
                  reduction = 'wnn.umap', max.cutoff = 2, 
                  cols = c("lightgrey","darkgreen"), ncol = 4)
DefaultAssay(ptall) <- 'RNA'
p7 <- FeaturePlot(ptall, features = c("CD3E", "CD8A", "PDCD1", "CTLA4"), 
                  reduction = 'wnn.umap', max.cutoff = 2, ncol = 4)
png("Featureplot_allcells_adt.png", width=40, height=15, units="cm", res=500)
p6
dev.off()
png("Featureplot_allcells_RNA.png", width=40, height=15, units="cm", res=500)
p7
dev.off()

#Find all markers
ptall.adtmarkers <- FindAllMarkers(ptall, min.pct = 0.15, logfc.threshold = 0.15, assay = "ADT")
ptall.rnamarkers <- FindAllMarkers(ptall, min.pct = 0.15, logfc.threshold = 0.15, assay = "RNA")

#Find top 20 markers
#adt
top20adt <- ptall.adtmarkers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
#rna
top20rna <- ptall.rnamarkers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
