library(Seurat)
library(dplyr)
library(scran)
library(harmony)
library(ggplot2)
library(cowplot)
library(data.table)
library(tibble)
library(CellChat)
library(ComplexHeatmap)

combined_seurat <- readRDS(~/combined_seurat.rds")

# Extract the data of each group
PTU_group <- subset(combined_seurat, subset = group == "PTU")
Control_group <- subset(combined_seurat, subset = group == "Control")
P_group <- subset(combined_seurat, subset = group == "30P")
LT4_group <- subset(combined_seurat, subset = group == "LT4")

# Create a CellChat object
cellchat_ptu <- createCellChat(object = PTU_group, group.by = "cell_type")
cellchat_control <- createCellChat(object = Control_group, group.by = "cell_type")
cellchat_30p <- createCellChat(object = P_group, group.by = "cell_type")
cellchat_LT4 <- createCellChat(object = LT4_group, group.by = "cell_type")

CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)

cellchat_ptu@DB <- subsetDB(CellChatDB)
cellchat_control@DB<- subsetDB(CellChatDB)
cellchat_30p@DB<- subsetDB(CellChatDB)
cellchat_LT4@DB<- subsetDB(CellChatDB)

cellchat_ptu <- subsetData(cellchat_ptu)
cellchat_control <- subsetData(cellchat_control)
cellchat_30p <- subsetData(cellchat_30p)
cellchat_LT4 <- subsetData(cellchat_LT4)

cellchat_ptu <- identifyOverExpressedGenes(cellchat_ptu)
cellchat_ptu <- identifyOverExpressedInteractions(cellchat_ptu)

cellchat_control <- identifyOverExpressedGenes(cellchat_control)
cellchat_control <- identifyOverExpressedInteractions(cellchat_control)

cellchat_ptu <- identifyOverExpressedGenes(cellchat_ptu)
cellchat_ptu <- identifyOverExpressedInteractions(cellchat_ptu)

cellchat_30p <- identifyOverExpressedGenes(cellchat_30p)
cellchat_30p <- identifyOverExpressedInteractions(cellchat_30p)

cellchat_LT4 <- identifyOverExpressedGenes(cellchat_LT4)
cellchat_LT4 <- identifyOverExpressedInteractions(cellchat_LT4)

cellchat_ptu <- computeCommunProb(cellchat_ptu)
cellchat_ptu <- filterCommunication(cellchat_ptu, min.cells = 10)  # Set the minimum cell count threshold
cellchat_ptu <- computeCommunProbPathway(cellchat_ptu)
cellchat_ptu <- aggregateNet(cellchat_ptu)

cellchat_ptu_Size <- as.numeric(table(cellchat_ptu@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_ptu@net$count, vertex.weight = cellchat_ptu_Size, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_ptu@net$weight, vertex.weight = cellchat_ptu_Size, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Compute the network centrality scores
cellchat_ptu <- netAnalysis_computeCentrality(cellchat_ptu, slot.name = "netP")
sig_pattern_PTU <- netAnalysis_signalingRole_heatmap(cellchat_ptu, pattern = "all")

ptu_pathways.show.all <- cellchat_ptu@netP$pathways

cellchat_control <- computeCommunProb(cellchat_control)
cellchat_control <- filterCommunication(cellchat_control, min.cells = 10)  # Set the minimum cell count threshold
cellchat_control <- computeCommunProbPathway(cellchat_control)
cellchat_control <- aggregateNet(cellchat_control)


cellchat_control_Size <- as.numeric(table(cellchat_control@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_control@net$count, vertex.weight = cellchat_control_Size, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_control@net$weight, vertex.weight = cellchat_control_Size, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


# Compute the network centrality scores
cellchat_control <- netAnalysis_computeCentrality(cellchat_control, slot.name = "netP")
sig_pattern_control <- netAnalysis_signalingRole_heatmap(cellchat_control, pattern = "all")


cellchat_30p <- computeCommunProb(cellchat_30p)
cellchat_30p <- filterCommunication(cellchat_30p, min.cells = 10)  # Set the minimum cell count threshold
cellchat_30p <- computeCommunProbPathway(cellchat_30p)
cellchat_30p <- aggregateNet(cellchat_30p)


cellchat_30p_Size <- as.numeric(table(cellchat_30p@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_30p@net$count, vertex.weight = cellchat_30p_Size, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_30p@net$weight, vertex.weight = cellchat_30p_Size, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Compute the network centrality scores
cellchat_30p <- netAnalysis_computeCentrality(cellchat_30p, slot.name = "netP")
sig_pattern_30p <- netAnalysis_signalingRole_heatmap(cellchat_30p, pattern = "all")


cellchat_LT4 <- computeCommunProb(cellchat_LT4)
cellchat_LT4 <- filterCommunication(cellchat_LT4, min.cells = 10)  # Set the minimum cell count threshold
cellchat_LT4 <- computeCommunProbPathway(cellchat_LT4)
cellchat_LT4 <- aggregateNet(cellchat_LT4)


cellchat_LT4_Size <- as.numeric(table(cellchat_LT4@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_LT4@net$count, vertex.weight = cellchat_LT4_Size, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_LT4@net$weight, vertex.weight = cellchat_LT4_Size, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Compute the network centrality scores
cellchat_LT4 <- netAnalysis_computeCentrality(cellchat_LT4, slot.name = "netP")
sig_pattern_LT4 <- netAnalysis_signalingRole_heatmap(cellchat_LT4, pattern = "all")


########################################################################################
#cellchat <- readRDS("~/cellchat.rds")
object.list <- list(Control = cellchat_control,PTU = cellchat_ptu,P=cellchat_30p,LT4=cellchat_LT4)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = T,group = c(1:4),measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = T, group = c(1:4),measure = "weight")
gg1 + gg2

par(mfrow = c(1,3), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T,measure = "count",comparison = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T,measure = "count",comparison = c(2,3))

netVisual_diffInteraction(cellchat, weight.scale = T,measure = "weight",comparison = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T,measure = "weight",comparison = c(2,3))
netVisual_diffInteraction(cellchat, weight.scale = T,measure = "weight",comparison = c(2,4))


netVisual_diffInteraction(cellchat, weight.scale = T,measure = "count",targets.use = 7,label.edge = T,comparison = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T,measure = "count",sources.use = 7,label.edge = T,comparison = c(1,2))


netVisual_diffInteraction(cellchat, weight.scale = T,measure = "count",targets.use = 7,label.edge = T,comparison = c(2,3))
netVisual_diffInteraction(cellchat, weight.scale = T,measure = "count",sources.use = 7,label.edge = T,comparison = c(2,3))

netVisual_diffInteraction(cellchat, weight.scale = T,measure = "count",targets.use = 7,label.edge = T,comparison = c(2,4))
netVisual_diffInteraction(cellchat, weight.scale = T,measure = "count",sources.use = 7,label.edge = T,comparison = c(2,4))


# Use the `saveRDS` function to save objects.
saveRDS(cellchat_control, file = "cellchat_control.rds")
saveRDS(cellchat_ptu, file = "cellchat_ptu.rds")
saveRDS(cellchat_30p, file = "cellchat_30p.rds")
saveRDS(cellchat_LT4, file = "cellchat_LT4.rds")
saveRDS(cellchat, file = "cellchat.rds")


