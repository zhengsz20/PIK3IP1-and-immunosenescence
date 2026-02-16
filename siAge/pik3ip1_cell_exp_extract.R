library(Seurat)
library(dplyr)
library(data.table)
library(tidyr)
data <- readRDS("E:/Wang_NI_2025/scRNA-seqProcessedLabelledObject.rds")

# Obtain all samples and cell types
samples <- unique(data@meta.data$sampleName)
cell_types <- unique(data@meta.data$secondary_type)

# Create the result matrix (samples × cell types, default filled with 0)
result_matrix <- matrix(
  0, 
  nrow = length(samples),
  ncol = length(cell_types),
  dimnames = list(samples, cell_types)
)

# Calculate the sample-level average expression of PIK3IP1 in each cell type.
for (ct in cell_types) {
  # Extract the subset of the current cell type
  ct_cells <- subset(data, subset = secondary_type == ct)
  
  #If the subset has no cells, skip it.
  if (ncol(ct_cells) == 0) next 
  
  # Extract the expression value of PIK3IP1
  expr_data <- LayerData(ct_cells, assay = "RNA", layer = "data")["PIK3IP1", ]
  
  # Calculate the mean by grouping according to the samples.
  sample_means <- tapply(
    expr_data,
    ct_cells@meta.data$sampleName,
    mean,
    na.rm = TRUE
  )
  
  # Fill the result matrix
  result_matrix[names(sample_means), ct] <- sample_means
}

# Convert to a data frame and format it.
result_df <- as.data.frame(result_matrix) %>%
  tibble::rownames_to_column(var = "sample") 

head(result_df, 3)
fwrite(result_df,"pik3ip1_exp_cell_type.txt",quote = F,sep = "\t",row.names = F)
