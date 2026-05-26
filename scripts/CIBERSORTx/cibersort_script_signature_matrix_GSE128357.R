library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)

# Charger les donnes
############################
# Organoid 1
org1 <- Read10X_h5("C:/Users/l25024469/Documents/Déconvolution cellulaire/1_Data_index/GSE128357/GSM3672210_Day-50_NM-organoid_filtered_gene_bc_matrices_h5.h5")
seurat1 <- CreateSeuratObject(org1)
seurat1$sample <- "organoid1"

# Organoid 2
org2 <- Read10X_h5("C:/Users/l25024469/Documents/Déconvolution cellulaire/1_Data_index/GSE128357/GSM4276650_Day-50_NM-organoid_2_filtered_feature_bc_matrix.h5")
seurat2 <- CreateSeuratObject(org2)
seurat2$sample <- "organoid2"
############################

# Merge
############################

mat1 <- seurat1@assays$RNA$counts
mat2 <- seurat2@assays$RNA$counts

# Harmonisation gènes
common_genes <- intersect(rownames(mat1), rownames(mat2))
mat1 <- mat1[common_genes, , drop = FALSE]
mat2 <- mat2[common_genes, , drop = FALSE]

counts_mat <- cbind(mat1, mat2)

# metadata
meta1 <- seurat1@meta.data
meta2 <- seurat2@meta.data
meta <- rbind(meta1, meta2)

# Reconstruction Seurat propre
seurat_obj <- CreateSeuratObject(counts = counts_mat)
seurat_obj$sample <- meta$sample
############################

# Qualite
############################

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

seurat_obj <- subset(
  seurat_obj,
  nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 10
)
############################

# Normalisation + clustering
############################

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

DimPlot(seurat_obj, label = TRUE)

############################

# MARKERS
############################

markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)

# Top markers
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50)

############################

# Annotation des types cellulaires
############################

seurat_obj$cell_type <- plyr::mapvalues(
  seurat_obj$seurat_clusters,
  from = as.character(0:18),
  to = c(
    "Activated_cytotoxic_T_cells_CTLA4+",          
    "Inflammatory_T_cells_TNFSF15+",               
    "Cycling_cells_S_phase_G2M",                   
    "Stromal_fibroblast_like_LUM+",                
    "Differentiated_epithelial_SPRR_S100A7+",      
    "Glial_like_cells_GFAP+",                      
    "GABAergic_neurons_NPY_GABRB2+",               
    "Mature_neurons_NEUROD6_LHX2+",                
    "Neural_progenitors_ASCL1_NEUROG1+",           
    "Highly_proliferative_cells_HIST_high",        
    "Specialized_epithelial_NPHS2_KCNJ1+",         
    "Inflammatory_myeloid_MMP9+",                  
    "Tolerogenic_dendritic_cells_HLA-G+",          
    "Epithelial_macrophage_like_KRT_MARCO+",       
    "Muscle_cells_CKM_MYH+",                       
    "Activated_macrophages_CCL3_MMP9+",            
    "Metabolic_epithelial_CYP_ABCB11+",            
    "Endothelial_cells_PECAM1_CD93+",              
    "Secretory_epithelial_KRT7_SPINK1+"))

seurat_obj$cell_type <- factor(seurat_obj$cell_type)
table(seurat_obj$seurat_clusters)

# Check
table(seurat_obj$cell_type)

############################

# Signature matrix pour CIBERSORTx
############################

# Expression normalisée
expr_mat <- GetAssayData(seurat_obj, layer ="counts")

cell_types <- unique(seurat_obj$cell_type)

# Moyenne par type cellulaire
signature_list <- lapply(cell_types, function(ct) {
  cells <- colnames(seurat_obj)[seurat_obj$cell_type == ct]
  Matrix::rowMeans(expr_mat[, cells, drop = FALSE])
})

signature_mat <- do.call(cbind, signature_list)
colnames(signature_mat) <- cell_types

############################

# Output
############################

signature_df <- data.frame(
  GeneSymbol = rownames(signature_mat),
  signature_mat,
  check.names = FALSE)

write.table(
  signature_df,
  file = "CIBERSORTx_signature_matrix_GSE128357.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE)

############################


