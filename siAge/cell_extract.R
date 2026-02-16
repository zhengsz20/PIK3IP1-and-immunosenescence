library(Seurat)
library(dplyr)
library(data.table)
data <- readRDS("E:/Wang_NI_2025/scRNA-seqProcessedLabelledObject.rds")

# Define the list of cell types to be extracted.
cell_types <- c("CD4T", "CD8T", "NK", "B", "pDC", "mDC", "B_BCR_GNLY", "Plasma")

# Cyclically extract and save
for (ct in cell_types) {
  subset_obj <- subset(data, ident = ct)
  saveRDS(subset_obj, paste0(ct, ".rds"))
}

# Extract two subclasses of Monocytes.
monocytes_cd14 <- subset(data, ident = "Monocytes_CD14")
monocytes_cd16 <- subset(data, ident = "Monocytes_CD16")

# Merge into a new object (retaining all data)
merged_monocytes <- merge(monocytes_cd14, monocytes_cd16)

# Reset the "ident" to a uniform type.
Idents(merged_monocytes) <- "Monocytes"

# write
saveRDS(merged_monocytes, "Monocytes.rds")
