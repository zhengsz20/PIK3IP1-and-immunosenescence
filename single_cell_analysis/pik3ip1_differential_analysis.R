library(Seurat)
library(dplyr)
library(tidyr)
combined_seurat <- readRDS("~/combined_seurat.rds")

# Extract metadata and add gene expression values
meta_data <- combined_seurat[[]] %>%
  mutate(Pik3ip1 = GetAssayData(combined_seurat, assay = "RNA", layer = "data")["Pik3ip1", ])

# Define the comparison function
perform_comparison <- function(ct, group1, group2) {
  subset_data <- meta_data %>%
    filter(cell_type == ct, group %in% c(group1, group2))
  
  # Calculate the average expression
  means <- subset_data %>%
    group_by(group) %>%
    summarise(mean_exp = mean(expm1(Pik3ip1))) %>%
    pivot_wider(names_from = group, values_from = mean_exp)
  
  # Wilcoxon rank sum test
  p_value <- wilcox.test(
    Pik3ip1 ~ group, 
    data = subset_data,
    alternative = "two.sided"
  )$p.value
  
  # Calculate Fold Change
  fc <- means[[group2]] / means[[group1]]
  
  # return to the result
  data.frame(
    cell_type = ct,
    comparison = paste0(group2, "_vs_", group1),
    mean_ref = means[[group1]],
    mean_target = means[[group2]],
    FC = fc,
    p_value = p_value
  )
}

#Perform all comparisons
results <- data.frame()
cell_types <- unique(meta_data$cell_type)

for (ct in cell_types) {
  # PTU vs Control
  if (sum(meta_data$cell_type == ct & meta_data$group == "Control") > 3 &&
      sum(meta_data$cell_type == ct & meta_data$group == "PTU") > 3) {
    results <- rbind(results, perform_comparison(ct, "Control", "PTU"))
  }
  
  # 30P vs PTU
  if (sum(meta_data$cell_type == ct & meta_data$group == "PTU") > 3 &&
      sum(meta_data$cell_type == ct & meta_data$group == "30P") > 3) {
    results <- rbind(results, perform_comparison(ct, "PTU", "30P"))
  }
  
  # LT4 vs PTU
  if (sum(meta_data$cell_type == ct & meta_data$group == "PTU") > 3 &&
      sum(meta_data$cell_type == ct & meta_data$group == "LT4") > 3) {
    results <- rbind(results, perform_comparison(ct, "PTU", "LT4"))
  }
}

# printf
final_table <- results %>%
  mutate(
    FC = signif(FC, 3),
    p_value = signif(p_value, 3),
    mean_ref = signif(mean_ref, 4),
    mean_target = signif(mean_target, 4)
  ) %>%
  arrange(cell_type, comparison)

# Check the results.
View(final_table)
final_table$adj_p <- p.adjust(final_table$p_value, method = "BH")
