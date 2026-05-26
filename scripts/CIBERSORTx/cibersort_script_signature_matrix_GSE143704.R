#ce script peut être utilisé uniquement dans le cas où au moment de trouver l'index su NCBI GEO, on peut obtenir le RAWDATA et le METADATA complets

library(data.table)

#Chargement data
########################
expr <- fread("C:/Users/l25024469/Documents/Deconvolution_cellulaire/1_Data_index/GSE143704_DeMicheli_HumanMuscleAtlas_rawdata.txt.gz",data.table=FALSE)
meta <- fread("C:/Users/l25024469/Documents/Deconvolution_cellulaire/1_Data_index/GSE143704_DeMicheli_HumanMuscleAtlas_metadata.txt.gz",data.table=FALSE)
########################

#Preparation des donnes
########################
colnames(expr)[1] <- "Gene"
rownames(expr) <- expr$Gene
expr$Gene <- NULL
########################

#Obtention des donnes de metadata
########################
cell_id <- meta$V1
cell_type <- meta$cell_annotation
########################

#Recuperer uniquement les genes qui sont presents dans les 2 datas
########################
keep <- cell_id %in% colnames(expr)
cell_id <- cell_id[keep]
cell_type <- cell_type[keep]
expr <- expr[, cell_id, drop = FALSE]
########################

#Construction de la signature matrix (on fait la moyenne par type cellulaire pour eviter repetition de IDs)
########################
signature <- sapply(unique(cell_type), function(ct) {
  rowMeans(expr[, cell_id[cell_type == ct], drop = FALSE])
})
signature <- as.data.frame(signature)
########################


#Export pour CIBERSORTx
########################
sig <- cbind(GeneSymbol = rownames(signature), signature)
rownames(sig) <- NULL

write.table(sig,
            "signature_matrix_GSE143704_CIBERSORTx_ready.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
########################