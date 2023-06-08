
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

FitGaussian <- function(norm_protein_expr) {
  npe <- norm_protein_expr[norm_protein_expr != 0]
  fit = Mclust(npe, G=2, model="V", verbose=FALSE)
  if (!is.null(fit)){
    signal <- as.numeric(which.max(fit$parameters$mean))
    expr_clean <- pnorm(norm_protein_expr,
                        mean = fit$parameters$mean[signal],
                        sd = sqrt(fit$parameters$variance$sigmasq[signal]))
  }else{
    
    expr_clean <- 0
  }
  
  return(expr_clean)
}


GaussCITE <-function (norm_cite_protein, num_cores = 1) 
{
  if (num_cores > 1) {
    cite_protein_list <- mclapply(1:ncol(norm_cite_protein), 
                                  function(i)  FitGaussian(norm_cite_protein[, i]), 
                                  mc.cores = num_cores)
  }
  else {
    cite_protein_list <- lapply(1:ncol(norm_cite_protein), 
                                function(i) FitGaussian(norm_cite_protein[, i]))
  }
  cite_protein_clean <- norm_cite_protein
  for (i in 1:ncol(norm_cite_protein)) {
    #      if (!is.null(cite_protein_list[[i]])){
    cite_protein_clean[, i] <- cite_protein_list[[i]]
    #     }
    #    else{
    #      print(i)
    #      cite_protein_clean[, i] <- 0
  }
  #}
  return(cite_protein_clean)
}


NormCells <- function(to_norm, norm_by=to_norm) {
  nonzero <- rowSums(norm_by) != 0
  protein_norm <- to_norm
  protein_norm[nonzero,] <- protein_norm[nonzero,] / rowSums(norm_by[nonzero,])
  return(protein_norm)
}

gaussNorm <- function (codex_filtered) 
{
  codex_clean <- codex_filtered
  for (i in 1:ncol(codex_filtered)) {
    fit = Mclust(codex_filtered[, i], G = 2, model = "V", 
                 verbose = FALSE)
    if (!is.null(fit)){
      signal <- as.numeric(which.max(fit$parameters$mean))
      #print(fit, fit$parameters$mean)
      expr_clean <- pnorm(codex_filtered[, i], mean = fit$parameters$mean[signal], 
                          sd = sqrt(fit$parameters$variance$sigmasq[signal]))
      codex_clean[, i] <- expr_clean
    }
    else{
      print(paste0("failed: ",i))
      codex_clean[, i] <- 0
    }
  }
  return(codex_clean)
}


option_list <- list (
  make_option(c("-i","--inPrefix"), type = 'character',
              help= "Path to folder containing codex data: _Prot.csv, and _Spat.csv files")
)

parser <- OptionParser(option_list = option_list)
arguments <- parse_args(parser, positional_arguments = FALSE)
print(arguments)

codex1_protein <- read.table(paste0(arguments$inPrefix,'_Prot.csv'), sep = ',', row.names =1,header=TRUE, stringsAsFactors = FALSE)
codex1_spatial <- read.table(paste0(arguments$inPrefix,'_Spat.csv'), sep = ',', row.names =1,header=TRUE, stringsAsFactors = FALSE)

codex1_size <- codex1_spatial$size
codex1_spatial_nm <- as.data.frame(cbind(x=codex1_spatial$x, y=codex1_spatial$y, z=codex1_spatial$z), row.names=row.names(codex1_spatial))
codex1_subset <- rownames(codex1_protein)[sample(nrow(codex1_protein),nrow(codex1_protein))]
codex1_proteins <- codex1_protein[codex1_subset,]
codex1_proteins <- codex1_proteins[, colSums(codex1_proteins) > 0]
codex1_blanks <- codex1_proteins[codex1_subset,FALSE]
codex1_blanks$BLANK0 <- 0
codex1_sizes <- codex1_spatial[codex1_subset,'size']
codex1_spatial_nms <- codex1_spatial_nm[codex1_subset,]
sizeMinMax <- c(min(codex_size)-1, max(codex_size)+1)

stvea1_objectRAW <- SetDataCODEX(codex_protein = codex1_proteins,
                                 codex_blanks = codex1_blanks,
                                 codex_size = codex1_sizes,
                                 codex_spatial = codex1_spatial_nms )

stvea1_objectRAW <- FilterCODEX(stvea1_objectRAW, size_lim = c(sizeMinMax[[1]],sizeMinMax[[2]]),
                                blank_lower = c(0, 0, 0, 0, 0, 0, 0),
                                blank_upper = c(15000,15000,15000,15000,15000,15000,15000))

stvea1_objectRAW@codex_emb <- NULL
stvea1_objectRAW@codex_clusters <- c('')
stvea1_objectRAW@codex_knn <- NULL
stvea1_objectRAW@codex_mRNA <- NULL
stvea1_objectRAW@corrected_codex <- NULL


codex_protein_norm <- stvea1_objectRAW@codex_protein - min(stvea1_objectRAW@codex_protein)
avg_cell_total <- mean(rowSums(codex_protein_norm))
codex_protein_norm <- NormCells(codex_protein_norm) * avg_cell_total
codex_clean <- gaussNorm(codex_protein_norm)
stvea1_objectRAW@codex_clean <-  codex_clean[, colSums(codex_clean) > 0]
stvea1_objectRAW@codex_protein <-  stvea1_objectRAW@codex_protein[, colSums(codex_clean) > 0]
stvea1_object <- copy(stvea1_objectRAW)



# CITE-seq protein
cite_protein <- read.table("rawdata/MIA_iADT.csv",
                           sep=",", row.names=1, header=TRUE)#, stringsAsFactors = FALSE)


cite_protein <- as.data.frame(cite_protein)
print(head(cite_protein))
cite_protein <- cite_protein[, !grepl("Mouse", colnames(cite_protein))]
cite_protein <- cite_protein[, !grepl("Rat", colnames(cite_protein))]
colnames(cite_protein) <- sapply(colnames(cite_protein), function(x) gsub("CITE.", "", x))
colnames(cite_protein)[names(cite_protein) == 'CD11C'] <- 'CD11c'
colnames(cite_protein)[names(cite_protein) == 'PDL1'] <- 'PD.L1'
#colnames(cite_protein)[names(cite_protein) == 'PD1'] <- 'PD.1'
#colnames(cite_protein)[names(cite_protein) == 'TIM3'] <- 'TIM.3'
colnames(cite_protein)[names(cite_protein) == 'CD8A'] <- 'CD8'
colnames(cite_protein)[names(cite_protein) == 'CD3'] <- 'CD3e'

#colnames(cite_protein)[names(cite_protein) == 'HLA.A.B.C'] <- 'HLA.X'
commonProt <- sort(intersect(colnames(cite_protein), colnames(stvea1_object@codex_clean)))
print(commonProt)
#sort(setdiff(colnames(cite_protein), colnames(colnames(stvea1_object@codex_clean))))
#sort(setdiff(colnames(stvea1_object@codex_clean), colnames(cite_protein)))
#cite_mRNA <- readRDS("rawdata/MIA_iRNA.RDS")
print("Reading RNA...")
cite_mRNA <- readRDS("rawdata/MIA_iRNA.RDS")

# CITE-seq mRNA latent space, such as output by scVI
cite_latent <- (read.table("rawdata/MIA_iPC_latent.csv",
                           sep=",", header=TRUE, row.names=1, stringsAsFactors = F)),                           
print("Done reading...")

cite_latent <- as.data.frame(cite_latent)
stvea1_object <- SetDataCITE(cite_mRNA= cite_mRNA,
                             cite_mRNA_norm = cite_mRNA,
                             cite_protein = cite_protein,
                             cite_latent = cite_latent,
                             stvea_object = stvea1_object)

print("Cleaning CITE-ADT...")
print(commonProt)
keepProt <- commonProt#[colSums(stvea1_object@cite_protein[,commonProt])>0]
print(keepProt)#write.csv(stvea1_object@cite_clean,"datasets/MIA_CiteClean.csv")
stvea1_object@cite_clean <- GaussCITE(stvea1_object@cite_protein[,keepProt])

print("Done cleaning cite-protein...")

stvea1_object@cite_mRNA_norm <- cite_mRNA
print(colSums(stvea1_object@cite_clean[,keepProt]))
stvea1_object@cite_clean <- stvea1_object@cite_clean[,keepProt]
stvea1_object@cite_protein <- stvea1_object@cite_protein[,keepProt]
print(keepProt)
print("Saving output...")

saveRDS(stvea1_object, paste0(arguments$inPrefix,'_Stvea.RDS'))
