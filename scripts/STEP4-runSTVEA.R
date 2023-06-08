
library(data.table)
library(mltools)
library(tidyr)
library(Matrix)
library(RColorBrewer) 
suppressPackageStartupMessages(library (mclust, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library (umap, warn.conflicts = FALSE, quietly = TRUE))

suppressPackageStartupMessages(library (Seurat, warn.conflicts = FALSE, quietly = TRUE))
suppressPackageStartupMessages(library (methods, warn.conflicts = FALSE, quietly = TRUE))
require(STvEA)
set.seed(20343)
print(.libPaths()) 
print(.packages())
suppressPackageStartupMessages(library (optparse, warn.conflicts = FALSE, quietly = TRUE))


option_list <- list (
    make_option(c("-i","--inPrefix"), type = 'character',
              help= "Path to folder containing codex data: _Prot.csv, and _Spat.csv files"),
    
    make_option(c("-n","--numChunks"), default = 10, type = 'integer',
                help= "Number of chunks to divide the query matrix into. default = %default .")
              )

parser <- OptionParser(option_list = option_list)
arguments <- parse_args(parser, positional_arguments = FALSE)
print(arguments)
nC <- as.integer(arguments$numChunks)
stvea1_object<- readRDS(paste0(arguments$inPrefix,'_Stvea.RDS'))



codex_split <- sample.int(nC,nrow(stvea1_object@codex_clean),replace=TRUE)

common_proteins <- colnames(stvea1_object@cite_clean)[colnames(stvea1_object@cite_clean) %in% colnames(stvea1_object@codex_clean)]
print(common_proteins)
ref_mat <- stvea1_object@cite_clean[,common_proteins]
stvea1_object@corrected_codex <- data.frame(matrix(,ncol=length(common_proteins),nrow=nrow(stvea1_object@codex_clean)))
rownames(stvea1_object@corrected_codex) <- row.names(stvea1_object@codex_clean)
colnames(stvea1_object@corrected_codex) <- common_proteins 
print(dim(stvea1_object@corrected_codex))
print(head(stvea1_object@corrected_codex))
for (ix in c(1:nC)){
  query_mat = stvea1_object@codex_clean[codex_split==ix,common_proteins]
  print(dim(query_mat))
  rna_mat = stvea1_object@cite_latent
  cite_index = 1
  num.cc = ncol(ref_mat)-1
  k.anchor = 20
  k.filter=100
  k.score=80
  k.weight=100
  verbose=TRUE
  
  # Call the same functions MapCODEXtoCITE() calls
  cca_matrix <- STvEA::RunCCA(t(ref_mat), t(query_mat), standardize=TRUE, num.cc=num.cc)$ccv
  neighbors <- FindNNrna(ref_emb = cca_matrix[1:nrow(ref_mat),],
                         query_emb = cca_matrix[(nrow(ref_mat)+1):nrow(cca_matrix),],
                         rna_mat = rna_mat,
                         cite_index = cite_index,
                         k=max(k.anchor, k.score), verbose=verbose)
  anchors <- FindAnchorPairs(neighbors, k.anchor=k.anchor)
  filteredAnchors <- FilterAnchors(ref_mat, query_mat, anchors, k.filter=k.filter, verbose=verbose)
  scoredAnchors <- ScoreAnchors(neighbors, filteredAnchors, nrow(ref_mat), nrow(query_mat), k.score=k.score, verbose=verbose)
  integration.matrix <- FindIntegrationMatrix(ref_mat, query_mat, neighbors, scoredAnchors, verbose=verbose)
  weights <- FindWeights(neighbors, scoredAnchors, query_mat, integration.matrix, k.weight=k.weight, verbose=verbose)
  corrected_data <- TransformDataMatrix(ref_mat, query_mat, integration.matrix, weights, verbose=verbose)
  
  # corrected_data is an rbind of ref_mat (cite_clean) and corrected_codex
  stSplit <- strsplit(arguments$inPrefix,"/")[[1]]
  saveRDS(corrected_data,  paste0(stSplit[1],'/',stSplit[2],'/','chunk',ix,'.RDS'))
  cN <- colnames(corrected_data)
  stvea1_object@corrected_codex[codex_split == ix,cN] <- corrected_data[(nrow(ref_mat)+1):nrow(corrected_data),]
  print(head(stvea1_object@corrected_codex[codex_split == ix,cN]))
}
saveRDS(stvea1_object, paste0(arguments$inPrefix,'_Output.RDS'))

saveRDS(codex_split,  paste0(arguments$inPrefix,'_split.RDS'))
