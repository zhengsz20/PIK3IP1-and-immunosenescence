library(magrittr)
library(oligo)
library(primeviewcdf)
library(pdInfoBuilder)
library(affy)
library(dplyr)

# Read the CEL file
raw_data <- ReadAffy(celfile.path = "F:/Framingham Heart Study_2013_15/GSE123696_RAW/")

# RMA normalization
eSet <- rma(raw_data)
expr_matrix <- exprs(eSet)  
head(expr_matrix[, 1:4])    


# Read the GPL15207 annotation file.
gpl_file <- "F:/Framingham Heart Study_2013_15/GPL15207-17536.txt"
annot_data <- read.delim(gpl_file, comment.char = "#", stringsAsFactors = FALSE, check.names = FALSE)

# Extract probe IDs and map gene symbols
probe_symbol <- annot_data %>%
  select(ID = `ID`, 
         symbol = `Gene Symbol`) %>%  # 列名根据实际文件调整
  filter(!is.na(symbol) & symbol != "")

# Merge the expression matrix with the gene symbols.
expr_annot <- expr_matrix %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("probe_id") %>%
  left_join(probe_symbol, by = c("probe_id" = "ID")) %>%
  filter(!is.na(symbol))  # 移除无基因符号的探针

# Handling duplicate gene symbols (taking the average of the expression values
expr_final <- expr_annot %>%
  group_by(symbol) %>%
  summarise(across(where(is.numeric), mean)) %>%
  tibble::column_to_rownames("symbol")

# View the expression matrix after annotation
head(expr_final[, 1:4])

# Convert the matrix to a data frame.
expr_df <- as.data.frame(expr_final)

# Add a gene column (extract row names and insert into the first column)
expr_df <- cbind(Gene = rownames(expr_df), expr_df)

# Delete the original row name
rownames(expr_df) <- NULL
head(expr_df[, 1:3]) 

fwrite(expr_df,"log2_GSE123698.txt",quote = F,sep = "\t",row.names = F)

# Extract the gene expression matrix from expr_df (excluding the Gene column).
expr_matrix <- as.matrix(expr_df[, -1])
rownames(expr_matrix) <- expr_df$Gene

# Calculate the Z statistic (standardized by gene)
z_matrix <- t(scale(t(expr_matrix)))

#Handling the case where the standard deviation is 0 (to avoid NaN)
z_matrix[is.nan(z_matrix)] <- 0

#Convert to a data frame and add a "Gene" column.
z_df <- as.data.frame(z_matrix)
z_df$Gene <- rownames(z_matrix)
rownames(z_df) <- NULL

# Adjust the column order (with the "Gene" column as the first column)
expr_z_df <- z_df[, c("Gene", colnames(z_matrix))]
fwrite(expr_z_df,"Z_GSE123698.txt",quote = F,sep = "\t",row.names = F)
