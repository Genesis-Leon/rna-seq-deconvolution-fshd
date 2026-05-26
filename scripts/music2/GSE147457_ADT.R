library(data.table) 
library(Matrix)
library(Seurat) 
library(tidyverse)
library(SingleCellExperiment)
library(MuSiC)

# Charger les donnes
############################
base_dir <- "~/GSE147457_RAW_ADT/"
files <- list.files(
  path = base_dir,
  pattern = "collapsed.HUMAN.dge.filtered.tsv.gz$",
  full.names = TRUE,
  recursive = TRUE)

length(files) #Only a verification
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

#Integration 
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
DimPlot(seurat_obj_a, label = TRUE)
# Vérification batch effect (samples)
DimPlot(seurat_obj_a, group.by = "sample")

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

# SINGLE CELL EXPERIMENT 
############################
counts_mat <- seurat_obj_a@assays$RNA$counts

sce <- SingleCellExperiment(
  assays = list(counts = counts_mat)
)

sce$cell_type <- seurat_obj_a$cell_type
sce$sample <- seurat_obj_a$sample

# CLEAN NA
sce <- sce[, !is.na(sce$cell_type)]
############################

# BULK DATA
############################
bulk <- fread("FSHD2_day30.txt")
bulk_clean <- bulk %>%
  group_by(Genes) %>%
  summarise(across(where(is.numeric), mean))

bulk_mat <- as.matrix(bulk_clean[, -1])
rownames(bulk_mat) <- bulk_clean$Genes

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
  file = "music_FSHD2_day30_GSE147457_ADT.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)



