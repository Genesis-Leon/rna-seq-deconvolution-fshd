library(Seurat)
library(SingleCellExperiment)
library(MuSiC)
library(dplyr)
library(data.table)

# Chargement de donnes
############################
# Organoid 1
org1 <- Read10X_h5("C:/Users/l25024469/Documents/Deconvolution_cellulaire/1_Data_index/GSE128357/GSM3672210_Day-50_NM-organoid_filtered_gene_bc_matrices_h5.h5")
seurat1 <- CreateSeuratObject(org1)
seurat1$sample <- "organoid1"

# Organoid 2
org2 <- Read10X_h5("C:/Users/l25024469/Documents/Deconvolution_cellulaire/1_Data_index/GSE128357/GSM4276650_Day-50_NM-organoid_2_filtered_feature_bc_matrix.h5")
seurat2 <- CreateSeuratObject(org2)
seurat2$sample <- "organoid2"
############################

# Merge
############################
mat1 <- seurat1@assays$RNA$counts # IMPORTANT: NE PAS TOUCHER avoid Seurat v5 multi-layer issues
mat2 <- seurat2@assays$RNA$counts

common_genes <- intersect(rownames(mat1), rownames(mat2))
mat1 <- mat1[common_genes, , drop = FALSE]
mat2 <- mat2[common_genes, , drop = FALSE]

counts_mat <- cbind(mat1, mat2)

meta <- rbind(seurat1@meta.data,seurat2@meta.data)

seurat_obj <- CreateSeuratObject(counts_mat)
seurat_obj$sample <- meta$sample
############################

# Qualite
############################
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

seurat_obj <- subset(
  seurat_obj,
  nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 10)
############################

# Normalisation + Clusterisation (Tiujours  faire)
############################
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
############################

# Annotation des types des cellules 
############################

markers <- FindAllMarkers( seurat_obj, only.pos = TRUE)  #Tous les marqueurs trouves dans chaque cluster
#Cette partie prend environ 1h donc il y a le temps de faire des autres choses :D

top_markers <- markers %>% #Top marqueurs utilises pour faire la annotaion
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

top_markers %>% arrange(cluster, desc(avg_log2FC)) #Voir les vraies marqeurs 

top_markers %>%
  group_by(cluster) %>%
  summarise(genes = paste(gene, collapse = ", "))

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

DimPlot(seurat_obj, label = TRUE, group.by = "seurat_clusters") #Visualisaton propre
############################

# SINGLE CELL EXPERIMENT 
############################

counts_mat <- seurat_obj@assays$RNA$counts #ATTENTION NE PAS CHANGER SINON PROBLEMES AVEC SEURAT V5

sce <- SingleCellExperiment(
  assays = list(counts = counts_mat)
)

sce$cell_type <- seurat_obj$cell_type
sce$sample <- seurat_obj$sample

# CLEAN NA
sce <- sce[, !is.na(sce$cell_type)]

############################

# BULK DATA
############################

bulk <- fread("C:/Users/l25024469/Documents/Deconvolution_cellulaire/1_Data/FSHD2/FSHD2_day30.txt")

bulk_clean <- bulk %>%
  group_by(gene) %>%
  summarise(across(where(is.numeric), mean))

bulk_mat <- as.matrix(bulk_clean[, -1])
rownames(bulk_mat) <- bulk_clean$gene

bulk_mat <- log2(bulk_mat + 1)

############################

# GENE INTERSECTION
############################

common_genes <- intersect(rownames(bulk_mat), rownames(counts_mat))

bulk_mat <- bulk_mat[common_genes, ]
counts_mat <- counts_mat[common_genes, ]

sce_1 <- sce[common_genes, ]

############################

# MUSIC
############################

result_music <- music_prop(
  bulk.mtx = bulk_mat,
  sc.sce = sce_1,
  clusters = "cell_type",
  samples = "sample"
)

############################

#OUTPUT
############################

prop <- result_music$Est.prop.weighted

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
  file = "music_FSHD2_day30_GSE128357.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


