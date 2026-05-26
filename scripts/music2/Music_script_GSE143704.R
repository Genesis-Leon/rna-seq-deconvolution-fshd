library(MuSiC)
library(SingleCellExperiment)
library(dplyr)
library(data.table)


#Bulk RNA-seq
#############
bulk <- fread("C:/Users/l25024469/Documents/Deconvolution_cellulaire/1_Data/FSHD2/FSHD2_day30.txt")

bulk_clean <- bulk %>% #Cette ligne permet de eviter un  decalage des resultats dans le output
  group_by(gene) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

bulk_mat <- as.matrix(bulk_clean[,-1])
rownames(bulk_mat) <- bulk_clean$gene

bulk_mat <- log2(bulk_mat + 1)
#############


#scRNA-seq
#############
sc_raw <- fread("C:/Users/l25024469/Documents/Deconvolution_cellulaire/1_Data_index/GSE143704_DeMicheli_HumanMuscleAtlas_rawdata.txt.gz",data.table = FALSE)

rownames(sc_raw) <- sc_raw[,1]
sc_raw <- as.matrix(sc_raw[,-1])

sc_meta <- fread("C:/Users/l25024469/Documents/Deconvolution_cellulaire/1_Data_index/GSE143704_DeMicheli_HumanMuscleAtlas_metadata.txt.gz",data.table = FALSE)

rownames(sc_meta) <- sc_meta[,1]
sc_meta <- sc_meta[,-1]
rownames(sc_meta) <- colnames(sc_raw)

stopifnot(all(rownames(sc_meta) == colnames(sc_raw)))
#############


#Gènes communs
#############
common_genes <- intersect(rownames(bulk_mat), rownames(sc_raw))

bulk_mat_filt <- bulk_mat[common_genes, ]
sc_raw_filt <- sc_raw[common_genes, ]
#############


#Créer SingleCellExperiment
#############
sc_sce <- SingleCellExperiment(
  assays = list(counts = sc_raw_filt),
  colData = sc_meta)

colData(sc_sce)$cell_annotation <- make.names(colData(sc_sce)$cell_annotation)
colData(sc_sce)$sampleID <- make.names(colData(sc_sce)$sampleID)
#############


#MuSiC
#############
res <- music_prop(
  bulk.mtx = bulk_mat_filt,
  sc.sce = sc_sce,
  clusters = "cell_annotation",
  samples = "sampleID",
  markers = NULL,
  normalize = TRUE)
#############


#Résultat
#############
prop <- res$Est.prop.weighted

if (is.list(prop)) {
  prop <- prop[[1]]
}

prop <- as.matrix(prop)

# forcer orientation correcte
if (nrow(prop) > ncol(prop)) {
  prop <- t(prop)
}

df_res17 <- as.data.frame(prop, check.names = FALSE)

df_res17$Sample <- rownames(df_res17)

df_res17 <- df_res17[, c("Sample", setdiff(colnames(df_res17), "Sample"))]

write.table(
  df_res17,
  file = "music_FSHD2_day30_GSE143704.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
#############

