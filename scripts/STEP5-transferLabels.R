require(STvEA)
library(data.table)
library(mltools)
library(tidyr)
library(Matrix)
suppressPackageStartupMessages(library (optparse, warn.conflicts = FALSE, quietly = TRUE))

#Rscript STEP5-transferLabels.R --inFile datasets/Wilmott_P7_Output.RDS --outPrefix datasets/Wilmott_P7
#-t rawdata/MIA_iUmapv0.csv -c leidenClust
GetTransferMatrixNew  <- function(from_dataset, to_dataset, k = floor(nrow(to_dataset)*0.01), c = 0.1) {
  # get corrected data for each dataset (either as input or from name in object)
  # compute query knn from CorNN
  # weight each nn based on gaussian kernel of distance
  # create weighted nn matrix as sparse matrix
  # return nn matrix (maybe as part of object)
  nn_list <- CorNN(to_dataset , from_dataset, k=k)
  nn_idx <- nn_list$nn.idx
  row.names(nn_idx) <- 1:nrow(nn_idx)
  
  nn_dists_exp <- exp(nn_list$nn.dists/-c)
  nn_weights <- nn_dists_exp / rowSums(nn_dists_exp)
  row.names(nn_weights) <- 1:nrow(nn_idx)
  
  sparse_coords <- gather(as.data.frame(t(nn_idx)))
  sparse_coords$key <- as.numeric(sparse_coords$key)
  sparse_entries <- gather(as.data.frame(t(nn_weights)))
  nn_matrix <- sparseMatrix(i = sparse_coords$value,
                            j = sparse_coords$key,
                            x = sparse_entries$value,
                            dims=c(nrow(to_dataset),nrow(from_dataset)))
  
  colnames(nn_matrix) <- row.names(from_dataset)
  row.names(nn_matrix) <- row.names(to_dataset)
  return(nn_matrix)
}

CorNN <- function(
  data,
  query = data,
  k = 10
) {
  t_data <- t(data)
  query <- as.matrix(query)
  neighbors <- matrix(rep(0, k*nrow(query)), ncol=k)
  distances <- matrix(rep(0, k*nrow(query)), ncol=k)
  for (i in 1:nrow(query)) {
    cor_dist <- 1-cor(query[i,], t_data)
    idx <- order(cor_dist)[1:k]
    neighbors[i,] <- idx
    distances[i,] <- cor_dist[idx]
  }
  return(list(nn.idx=neighbors, nn.dists=distances))
}

option_list <- list (
  make_option(c("-i","--inFile"), type = 'character',
              help= "Path to STvEA RDS object."),
  make_option(c("-t","--transferFile"), type = 'character',
              help= "Path to STvEA RDS object. Leave it empty to only transfer CITE-ADT and CITE-RNA."),
  make_option(c("-c","--columnName"), type = 'character',default = 'expression',
              help= "Column name in trasferFile to map. Leave it empty to only transfer CITE-ADT and CITE-RNA"),
  make_option(c("-o","--outPrefix"), type = 'character',
              help= "Prefix to output csv file.")
)

parser <- OptionParser(option_list = option_list)
arguments <- parse_args(parser, positional_arguments = FALSE)
print(arguments)
print("Reading input data.")
stvea1_object <- readRDS(arguments$inFile)
labelsDF <- read.csv(arguments$transferFile, sep=',', row.names = 1)
#print(head(labelsDF))
# to-do pass k and c as input arguments.
print("Computing transfer matrix. This may take a while.")
codex1_subset <- rownames(stvea1_object@corrected_codex)[sample(nrow(stvea1_object@corrected_codex),25000)]

transfer_matrix1 <- GetTransferMatrixNew(stvea1_object@corrected_codex[codex1_subset,], 
                                         stvea1_object@cite_clean[,colnames(stvea1_object@corrected_codex)])#,k = 15,  c=0.1)
saveRDS(transfer_matrix1, paste0(arguments$outPrefix,'_','transMatrix','.RDS'))

#transfer_matrix1
print("Done computing transfer matrix")
if (arguments$columnName != 'expression'){
# transfer specific column
cName <- arguments$columnName
labelsDF[,cName] <- as.factor(labelsDF[,cName])
#this will convert the cluster ids to cName_<clusterID>; which needs to be removed later.
onh <- one_hot(as.data.table(labelsDF),cols=cName)
#onh <- onh0[,grepl(cName, colnames(onh)),with = FALSE]

tfdf <- onh[,grepl(cName, colnames(onh)),with = FALSE]
#print(head(tfdf))
#print(dim(tfdf))
#print(dim(transfer_matrix1))
tflab <- as.matrix(t(transfer_matrix1))  %*%  Matrix(as.matrix(tfdf), sparse = T)
print(sum(is.na(tflab)))
# Retreive the cluster labels back by removing cName_ prefix.
TransferredClusters <- gsub(paste0(cName,"_"),"",colnames(tflab)[max.col(tflab,ties.method="first")])
#TransferredClusters <- colnames(tflab)[max.col(tflab,ties.method="first")]
outDF <- as.data.frame(stvea1_object@codex_clean[codex1_subset,])
# t_ represents transferred values.
outDF[,paste0("t_",cName)] <-TransferredClusters
row.names(outDF) <- row.names(stvea1_object@corrected_codex[codex1_subset,])
write.csv(outDF, paste0(arguments$outPrefix,'_',cName,'.csv'))
}
print("Done labels")
cADTNorm <- read.csv('../Swarbrick/Tsl_ADT_Raw_GMM.csv',sep=",", row.names=1, header=TRUE)
tfprot <- as.matrix(t(transfer_matrix1)  %*%  Matrix(as.matrix(cADTNorm), sparse = T))
protDF <- as.data.frame(tfprot)
#colnames(protDF) <- paste("t",colnames(stvea1_object@cite_clean),sep="_")
colnames(protDF) <- paste("t",colnames(cADTNorm),sep="_")
row.names(protDF) <- row.names(stvea1_object@corrected_codex[codex1_subset,])

write.csv(protDF, paste0(arguments$outPrefix,'_','cADT','.csv'))
# skipping this for now
print("Done ADT")

tfRNA <- as.matrix(t(transfer_matrix1)  %*%  Matrix(as.matrix(stvea1_object@cite_mRNA_norm), sparse = T))
rnaDF <- as.data.frame(tfRNA)
colnames(rnaDF) <- paste("t",colnames(stvea1_object@cite_mRNA_norm),sep="_")
row.names(rnaDF) <- row.names(stvea1_object@corrected_codex[codex1_subset,])
saveRDS(rnaDF, paste0(arguments$outPrefix,'_','cRNA','.RDS'))
write.csv(rnaDF, paste0(arguments$outPrefix,'_','cRNA','.csv'))

print("Done.")

