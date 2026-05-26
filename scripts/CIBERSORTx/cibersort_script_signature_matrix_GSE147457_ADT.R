library(data.table) 
library(Matrix) 
library(Seurat) 
library(tidyverse) 
library(SingleCellExperiment) 
library(MuSiC)

# Charger les donnes
############################

#Verifier où on a les données
base_dir <- "~/Bureau/datacenter/Eq_Magdinier/Etudiants/Genesis/Deconvolution_cellulaire/1_Data_index/GSE147457_RAW_ADT/" 
files <- list.files(
  path = base_dir,
  pattern = "collapsed.HUMAN.dge.filtered.tsv.gz$",
  full.names = TRUE,
  recursive = TRUE)

length(files) # change selon le dossier choisi 
############################

# Creation objet seurat
############################
seurat_list <- list()

for (f in files) {
  
  cat("Loading:", f, "\n")
  
  mat <- fread(f, data.table = FALSE)
  
  # Format classique GEO
  rownames(mat) <- mat[,1]
  mat <- mat[,-1]
  
  mat <- as.matrix(mat)
  mat <- Matrix(mat, sparse = TRUE)
  
  # Nom propre du sample
  sample_name <- sub(".collapsed.HUMAN.dge.filtered.tsv.gz", "", basename(f))
  
  seu <- CreateSeuratObject(counts = mat)
  seu$sample <- sample_name
  
  seurat_list[[sample_name]] <- seu
}
############################

#Integration PAS MERGE CAR IL Y A TROP DE FICHIERS
############################
# Normalisation 
seurat_list <- lapply(seurat_list, NormalizeData)
seurat_list <- lapply(seurat_list, FindVariableFeatures, nfeatures = 2000)

# PCA
seurat_list <- lapply(seurat_list, function(x) {
  x <- ScaleData(x, verbose = FALSE)
  x <- RunPCA(x, npcs = 30, verbose = FALSE)
  return(x)
})
#Features
features <- SelectIntegrationFeatures(seurat_list, nfeatures = 2000)
# Anchors
anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  anchor.features = features,
  reduction = "rpca",
  dims = 1:30,
  k.anchor = 5
)

# Integration
seurat_obj_a <- IntegrateData(anchorset = anchors)
############################

#clustering
############################
# Assay intégré pour clustering
DefaultAssay(seurat_obj_a) <- "integrated"
# Scaling
seurat_obj_a <- ScaleData(seurat_obj_a, verbose = FALSE)
# PCA
seurat_obj_a <- RunPCA(seurat_obj_a, npcs = 30, verbose = FALSE)
# Graph + clustering
seurat_obj_a <- FindNeighbors(seurat_obj_a, dims = 1:30)
seurat_obj_a <- FindClusters(seurat_obj_a, resolution = 0.5)
# UMAP
seurat_obj_a <- RunUMAP(seurat_obj_a, dims = 1:30)
# Visualisation clusters
DimPlot(seurat_obj_a, label = TRUE) #Les clusters sont bien separes donc on voit des groupes biologiquements bien separes 
# Vérification batch effect
DimPlot(seurat_obj_a, group.by = "sample") #Tous les samples sont bien melanges, donc il na pas de effet batch entre eux

DefaultAssay(seurat_obj_a) <- "RNA"
############################

# MARKERS
############################
seurat_obj_a <- JoinLayers(seurat_obj_a)
markers <- FindAllMarkers(seurat_obj_a, only.pos = TRUE)

# Top markers
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 20) 

top_markers %>% arrange(cluster, desc(avg_log2FC))

top_markers %>%
  group_by(cluster) %>%
  summarise(genes = paste(gene, collapse = ", "))

seurat_obj_a$cell_type <- plyr::mapvalues(
  seurat_obj_a$seurat_clusters,
  from = as.character(0:8),
  to = c(
    "Mesenchymal_stromal_cells_SPON1_HAND2",
    "Developmental_mesenchymal_progenitor_cells_HES7_WISP1",
    "Neuromuscular_cells_CDH15_CHRND",
    "Cycling_progenitor_like_cells_WNT5B_PCNA",
    "Fibroblast_like_cells_DKK1_FRAS1",
    "Metabolic_mesenchymal_cells_IBSP_CIDEA",
    "Endothelial_cells_RNASE1_PLVAP",
    "Peripheral_neuron_like_cells_SYT1_PRP",
    "Muscle_cells_MYOG_MYH6"
  )
)

seurat_obj_a$cell_type <- factor(seurat_obj_a$cell_type)
DimPlot(seurat_obj_a, label = TRUE, group.by = "seurat_clusters") 
############################

# Signature matrix pour CIBERSORTx
############################

# Expression normalisée
expr_mat <- GetAssayData(seurat_obj_a, layer ="counts")

cell_types <- unique(seurat_obj_a$cell_type)

# Moyenne par type cellulaire
signature_list <- lapply(cell_types, function(ct) {
  cells <- colnames(seurat_obj_a)[seurat_obj_a$cell_type == ct]
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
  file = "CIBERSORTx_signature_matrix_GSE147457_ADT.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE)

############################


