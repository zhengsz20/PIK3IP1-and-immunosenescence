# ---- load packages ----
library(Seurat)
library(dplyr)
library(scran)
library(harmony)
library(ggplot2)
library(cowplot)
library(data.table)
library(tibble)


# 
# Data loading

C30 <- Read10X(data.dir = "F:/单细胞/CellRanger_data/C30/")
C46 <- Read10X(data.dir = "F:/单细胞/CellRanger_data/C46/")
C65 <- Read10X(data.dir = "F:/单细胞/CellRanger_data/C65/")
C67 <- Read10X(data.dir = "F:/单细胞/CellRanger_data/C67/")


M61 <- Read10X(data.dir = "F:/单细胞/CellRanger_data/M61/")
M84 <- Read10X(data.dir = "F:/单细胞/CellRanger_data/M84/")
M88 <- Read10X(data.dir = "F:/单细胞/CellRanger_data/M88/")
M95 <- Read10X(data.dir = "F:/单细胞/CellRanger_data/M95/")

TP12 <- Read10X(data.dir = "F:/单细胞/CellRanger_data/TP12/")
TP16 <- Read10X(data.dir = "F:/单细胞/CellRanger_data/TP16/")
TP22 <- Read10X(data.dir = "F:/单细胞/CellRanger_data/TP22/")
TP31 <- Read10X(data.dir = "F:/单细胞/CellRanger_data/TP31/")

L9 <- Read10X(data.dir = "F:/单细胞/CellRanger_data/L9/")
L10 <- Read10X(data.dir = "F:/单细胞/CellRanger_data/L10/")
L53 <- Read10X(data.dir = "F:/单细胞/CellRanger_data/L53/")
L57<- Read10X(data.dir = "F:/单细胞/CellRanger_data/L57/")


# Creation of Seurat object

C30_seurat <- CreateSeuratObject(counts = C30, project = "C30")
C46_seurat <- CreateSeuratObject(counts = C46, project = "C46")
C65_seurat <- CreateSeuratObject(counts = C65, project = "C65")
C67_seurat <- CreateSeuratObject(counts = C67, project = "C67")

M61_seurat <- CreateSeuratObject(counts = M61, project = "M61")
M84_seurat <- CreateSeuratObject(counts = M84, project = "M84")
M88_seurat <- CreateSeuratObject(counts = M88, project = "M88")
M95_seurat <- CreateSeuratObject(counts = M95, project = "M95")

TP12_seurat <- CreateSeuratObject(counts = TP12, project = "TP12")
TP16_seurat <- CreateSeuratObject(counts = TP16, project = "TP16")
TP22_seurat <- CreateSeuratObject(counts = TP22, project = "TP22")
TP31_seurat <- CreateSeuratObject(counts = TP31, project = "TP31")


L9_seurat <- CreateSeuratObject(counts = L9, project = "L9")
L10_seurat <- CreateSeuratObject(counts = L10, project = "L10")
L53_seurat <- CreateSeuratObject(counts = L53, project = "L53")
L57_seurat <- CreateSeuratObject(counts = L57, project = "L57")

# free the memory
rm(C30,C46,C65,C67,M61, M84, M88, M95, TP12, TP16, TP22,TP31,L9,L10,L53,L57)
gc()

# Genes expressed in less than three cells were removed.

filter_features <- function(seurat_obj) {
  expressed_genes <- rownames(seurat_obj)[rowSums(GetAssayData(seurat_obj, layer = "counts") > 0) >= 3]
  subset(seurat_obj, features = expressed_genes)
}

C30_seurat <- filter_features(C30_seurat)
C46_seurat <- filter_features(C46_seurat)
C65_seurat <- filter_features(C65_seurat)
C67_seurat <- filter_features(C67_seurat)

M61_seurat <- filter_features(M61_seurat)
M84_seurat <- filter_features(M84_seurat)
M88_seurat <- filter_features(M88_seurat)
M95_seurat <- filter_features(M95_seurat)

TP12_seurat <- filter_features(TP12_seurat)
TP16_seurat <- filter_features(TP16_seurat)
TP22_seurat <- filter_features(TP22_seurat)
TP31_seurat <- filter_features(TP31_seurat)

L9_seurat <- filter_features(L9_seurat)
L10_seurat <- filter_features(L10_seurat)
L53_seurat <- filter_features(L53_seurat)
L57_seurat <- filter_features(L57_seurat)



# Calculate the proportion of mitochondrial genes and the proportion of red blood cells for each Seurat object.
C30_seurat[["percent.mt"]] <- PercentageFeatureSet(C30_seurat, pattern = "^MT-")
C46_seurat[["percent.mt"]] <- PercentageFeatureSet(C46_seurat, pattern = "^MT-")
C65_seurat[["percent.mt"]] <- PercentageFeatureSet(C65_seurat, pattern = "^MT-")
C67_seurat[["percent.mt"]] <- PercentageFeatureSet(C67_seurat, pattern = "^MT-")

M61_seurat[["percent.mt"]] <- PercentageFeatureSet(M61_seurat, pattern = "^MT-")
M84_seurat[["percent.mt"]] <- PercentageFeatureSet(M84_seurat, pattern = "^MT-")
M88_seurat[["percent.mt"]] <- PercentageFeatureSet(M88_seurat, pattern = "^MT-")
M95_seurat[["percent.mt"]] <- PercentageFeatureSet(M95_seurat, pattern = "^MT-")

TP12_seurat[["percent.mt"]] <- PercentageFeatureSet(TP12_seurat, pattern = "^MT-")
TP16_seurat[["percent.mt"]] <- PercentageFeatureSet(TP16_seurat, pattern = "^MT-")
TP22_seurat[["percent.mt"]] <- PercentageFeatureSet(TP22_seurat, pattern = "^MT-")
TP31_seurat[["percent.mt"]] <- PercentageFeatureSet(TP31_seurat, pattern = "^MT-")

L9_seurat[["percent.mt"]] <- PercentageFeatureSet(L9_seurat, pattern = "^MT-")
L10_seurat[["percent.mt"]] <- PercentageFeatureSet(L10_seurat, pattern = "^MT-")
L53_seurat[["percent.mt"]] <- PercentageFeatureSet(L53_seurat, pattern = "^MT-")
L57_seurat[["percent.mt"]] <- PercentageFeatureSet(L57_seurat, pattern = "^MT-")

C30_seurat[["percent.rbc"]] <- PercentageFeatureSet(C30_seurat, pattern = "^Hbb-")
C46_seurat[["percent.rbc"]] <- PercentageFeatureSet(C46_seurat, pattern = "^Hbb-")
C65_seurat[["percent.rbc"]] <- PercentageFeatureSet(C65_seurat, pattern = "^Hbb-")
C67_seurat[["percent.rbc"]] <- PercentageFeatureSet(C67_seurat, pattern = "^Hbb-")

M61_seurat[["percent.rbc"]] <- PercentageFeatureSet(M61_seurat, pattern = "^Hbb-")
M84_seurat[["percent.rbc"]] <- PercentageFeatureSet(M84_seurat, pattern = "^Hbb-")
M88_seurat[["percent.rbc"]] <- PercentageFeatureSet(M88_seurat, pattern = "^Hbb-")
M95_seurat[["percent.rbc"]] <- PercentageFeatureSet(M95_seurat, pattern = "^Hbb-")

TP12_seurat[["percent.rbc"]] <- PercentageFeatureSet(TP12_seurat, pattern = "^Hbb-")
TP16_seurat[["percent.rbc"]] <- PercentageFeatureSet(TP16_seurat, pattern = "^Hbb-")
TP22_seurat[["percent.rbc"]] <- PercentageFeatureSet(TP22_seurat, pattern = "^Hbb-")
TP31_seurat[["percent.rbc"]] <- PercentageFeatureSet(TP31_seurat, pattern = "^Hbb-")

L9_seurat[["percent.rbc"]] <- PercentageFeatureSet(L9_seurat, pattern = "^Hbb-")
L10_seurat[["percent.rbc"]] <- PercentageFeatureSet(L10_seurat, pattern = "^Hbb-")
L53_seurat[["percent.rbc"]] <- PercentageFeatureSet(L53_seurat, pattern = "^Hbb-")
L57_seurat[["percent.rbc"]] <- PercentageFeatureSet(L57_seurat, pattern = "^Hbb-")

C30_seurat[["percent.rbc2"]] <- PercentageFeatureSet(C30_seurat, pattern = "^Hba-")
C46_seurat[["percent.rbc2"]] <- PercentageFeatureSet(C46_seurat, pattern = "^Hba-")
C65_seurat[["percent.rbc2"]] <- PercentageFeatureSet(C65_seurat, pattern = "^Hba-")
C67_seurat[["percent.rbc2"]] <- PercentageFeatureSet(C67_seurat, pattern = "^Hba-")


M61_seurat[["percent.rbc2"]] <- PercentageFeatureSet(M61_seurat, pattern = "^Hba-")
M84_seurat[["percent.rbc2"]] <- PercentageFeatureSet(M84_seurat, pattern = "^Hba-")
M88_seurat[["percent.rbc2"]] <- PercentageFeatureSet(M88_seurat, pattern = "^Hba-")
M95_seurat[["percent.rbc2"]] <- PercentageFeatureSet(M95_seurat, pattern = "^Hba-")

TP12_seurat[["percent.rbc2"]] <- PercentageFeatureSet(TP12_seurat, pattern = "^Hba-")
TP16_seurat[["percent.rbc2"]] <- PercentageFeatureSet(TP16_seurat, pattern = "^Hba-")
TP22_seurat[["percent.rbc2"]] <- PercentageFeatureSet(TP22_seurat, pattern = "^Hba-")
TP31_seurat[["percent.rbc2"]] <- PercentageFeatureSet(TP31_seurat, pattern = "^Hba-")


L9_seurat[["percent.rbc2"]] <- PercentageFeatureSet(L9_seurat, pattern = "^Hba-")
L10_seurat[["percent.rbc2"]] <- PercentageFeatureSet(L10_seurat, pattern = "^Hba-")
L53_seurat[["percent.rbc2"]] <- PercentageFeatureSet(L53_seurat, pattern = "^Hba-")
L57_seurat[["percent.rbc2"]] <- PercentageFeatureSet(L57_seurat, pattern = "^Hba-")


# Quality control: Remove low-quality cells
filter_cells <- function(seurat_obj) {
  subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 2000 & nCount_RNA < 25000 & percent.mt < 10&percent.rbc<1&percent.rbc2<1)
}

C30_seurat <- filter_cells(C30_seurat)
C46_seurat <- filter_cells(C46_seurat)
C65_seurat <- filter_cells(C65_seurat)
C67_seurat <- filter_cells(C67_seurat)


M61_seurat <- filter_cells(M61_seurat)
M84_seurat <- filter_cells(M84_seurat)
M88_seurat <- filter_cells(M88_seurat)
M95_seurat <- filter_cells(M95_seurat)

TP12_seurat <- filter_cells(TP12_seurat)
TP16_seurat <- filter_cells(TP16_seurat)
TP22_seurat <- filter_cells(TP22_seurat)
TP31_seurat <- filter_cells(TP31_seurat)

L9_seurat <- filter_cells(L9_seurat)
L10_seurat <- filter_cells(L10_seurat)
L53_seurat <- filter_cells(L53_seurat)
L57_seurat <- filter_cells(L57_seurat)


# Add labels to each group.
C30_seurat$group <- "Control"
C46_seurat$group <- "Control"
C65_seurat$group <- "Control"
C67_seurat$group <- "Control"


M61_seurat$group <- "PTU"
M84_seurat$group <- "PTU"
M88_seurat$group <- "PTU"
M95_seurat$group <- "PTU"

TP12_seurat$group <- "30P"
TP16_seurat$group <- "30P"
TP22_seurat$group <- "30P"
TP31_seurat$group <- "30P"

L9_seurat$group <- "LT4"
L10_seurat$group <- "LT4"
L53_seurat$group <- "LT4"
L57_seurat$group <- "LT4"


# Merge the samples into one Seurat object
combined_seurat <- merge(
  C30_seurat,
  y = list(C46_seurat, C65_seurat, C67_seurat,M61_seurat,M84_seurat, M88_seurat, M95_seurat,TP12_seurat,TP16_seurat,TP22_seurat,TP31_seurat,L9_seurat,L10_seurat,L53_seurat,L57_seurat), 
  add.cell.ids = c("C30", "C46", "C65", "C67","M61", "M84", "M88", "M95","TP12", "TP16", "TP22", "TP31","L9","L10","L53","L57")
)
combined_seurat <- JoinLayers(object = combined_seurat,assay = "RNA",layers = "counts")

rm(C30_seurat,C46_seurat,C65_seurat,C67_seurat,M61_seurat,M84_seurat, M88_seurat, M95_seurat,TP12_seurat,TP16_seurat,TP22_seurat,TP31_seurat,L9_seurat,L10_seurat,L53_seurat,L57_seurat)

gc()


  # Draw the violin plot of nFeature_RNA
  p1 <- VlnPlot(seurat_obj, features = "nFeature_RNA", pt.size = 0, layer = "counts") + 
    labs(title = paste(title_prefix, "nFeature_RNA")) + 
    xlab("Sample") + 
    ylab("Number of Genes") + 
    theme_minimal_base
  
  # Draw the violin plot of nCount_RNA
  p2 <- VlnPlot(seurat_obj, features = "nCount_RNA", pt.size = 0, layer = "counts") + 
    labs(title = paste(title_prefix, "nCount_RNA")) + 
    xlab("Sample") + 
    ylab("UMI Counts") + 
    theme_minimal_base
  
  # Draw a violin plot for percent.mt
  p3 <- VlnPlot(seurat_obj, features = "percent.mt", pt.size = 0, layer = "counts") + 
    labs(title = paste(title_prefix, "percent.mt")) + 
    xlab("Sample") + 
    ylab("Percent Mitochondrial") + 
    theme_minimal_base
  
  # Merge three violin plots using the cowplot package.
  plot_combined <- plot_grid(p1, p2, p3, labels = c("A", "B", "C"), ncol = 1)
  
  return(plot_combined)
}

# Call the function to draw the quality control violin plot.
plot_combined <- plot_qc_violin(combined_seurat, title_prefix = "Combined QC")

# Display graphics
print(plot_combined)

# Data normalization
combined_seurat <- NormalizeData(combined_seurat)

# Identify highly variable genes
combined_seurat <- FindVariableFeatures(combined_seurat, selection.method = "vst", nfeatures = 4000)

# Perform PCA using high-variant genes
combined_seurat <- ScaleData(combined_seurat, features = VariableFeatures(combined_seurat))
combined_seurat <- RunPCA(combined_seurat, features = VariableFeatures(combined_seurat))

# Draw the elbow plot of PCA
ElbowPlot(combined_seurat,ndims=50)

# Harmony correction
combined_seurat <- RunHarmony(combined_seurat, group.by.vars = "orig.ident", dims = 1:10)

# Clustering and UMAP dimensionality reduction
combined_seurat <- RunUMAP(combined_seurat, reduction = "harmony", dims = 1:10)
combined_seurat <- FindNeighbors(combined_seurat, reduction = "harmony", dims = 1:10)

combined_seurat <- FindClusters(combined_seurat, resolution = 0.1)

# Obtain UMAP coordinates and clustering information
umap_data <- Embeddings(combined_seurat, reduction = "umap")
cluster_data <- combined_seurat$seurat_clusters

# Create a data frame containing UMAP coordinates and clustering information.
umap_df <- data.frame(
  UMAP_1 = umap_data[, 1],
  UMAP_2 = umap_data[, 2],
  Cluster = as.factor(cluster_data)
)

# Calculate the center point of each cluster.
cluster_centers <- umap_df %>%
  group_by(Cluster) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

# Draw a UMAP plot and label the cluster numbers.
umap_plot <- DimPlot(combined_seurat, reduction = "umap", group.by = "seurat_clusters", label = FALSE) +
  geom_text(data = cluster_centers, aes(x = UMAP_1, y = UMAP_2, label = Cluster), size = 5, color = "black")

# Display the UMAP plot
print(umap_plot)


# Cell cluster annotation: Manually assigned based on known marker genes.
# Create a named vector for annotation.
cluster_annotations <- c(
  "0" = "CD4 T cell",
  "1" = "CD8 T cell",
  "2" = "NK cell",
  "3" = "B cell",
  "4" = "nonclassical monocytes",
  "5" = "CD14 monocytes",
  "6" = "Endothelial cell"
)

# Add annotations to the Seurat object.
combined_seurat$cell_type <- factor(
  combined_seurat$seurat_clusters, 
  levels = names(cluster_annotations), 
  labels = cluster_annotations
)

# Visualized UMAP plot with annotations
DimPlot(combined_seurat, reduction = "umap", group.by = "cell_type", label = TRUE, label.size = 4)

#  Visualize the expression of a certain gene on UMAP
# eplace 'GENE_NAME' with the name of the gene you are interested in.
gene_of_interest <- "GENE_NAME"  
FeaturePlot(combined_seurat, features = gene_of_interest, reduction = "umap", cols = c("lightgrey", "blue"))




