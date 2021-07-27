library(dplyr)

############## ETV1 cell exploration #################

DimPlot(par.combined)
etv1.cells <- subset(par.combined, idents = "Etv1+")

DefaultAssay(etv1.cells) <- "RNA"
etv1.cells <- SplitObject(etv1.cells, split.by = "Treatment")

# NORMALIZE
etv1.cells <- lapply(X = etv1.cells, FUN = function(y){
  y <- SCTransform(y, vars.to.regress = "percent.mt", verbose = T)
})

features <- SelectIntegrationFeatures(etv1.cells, nfeatures = 3000)
etv1.cells <- PrepSCTIntegration(etv1.cells, anchor.features = features)
anchors <- FindIntegrationAnchors(etv1.cells, anchor.features = features, normalization.method = "SCT")
etv1.cells <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = 70)
etv1.cells <- RunPCA(etv1.cells, npcs = 30)
etv1.cells <- RunUMAP(etv1.cells, reduction = "pca", dims = 1:5)
etv1.cells <- FindNeighbors(etv1.cells, dims = 1:5, force.recalc = T)
etv1.cells <- FindClusters(etv1.cells, resolution = 0.4 )
DimPlot(etv1.cells, group.by = "Treatment")

etv1.sub.markers <- FindAllMarkers(etv1.cells, only.pos = T)
etv1.sub.markers <- etv1.sub.markers %>% filter(p_val_adj <0.05)
etv1.sub.markers.top <- etv1.sub.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

FeaturePlot(etv1.cells, features = "Ccnd1",
            min.cutoff = "q10"
            )

etv1matrix <- AverageExpression(etv1.cells, features = unique(etv1.sub.markers.top$gene), assays = "SCT")
etv1matrix <- etv1matrix$SCT

pheatmap::pheatmap(etv1matrix, scale = "row", fontsize = 12, cluster_rows = F, cluster_cols = F,
                   border_color = "black")


# ============================================================ #
# integrate with SMG gland to see what these cells are similar to #


SMG <- readRDS("D:/My GitHub/scRNAseq_mousePG/Data/Postnatal SMG Integrated.rds")

DefaultAssay(SMG) <- "RNA"
SMG <- SplitObject(SMG, split.by = "stage")


parotid <- par.combined

DefaultAssay(parotid) <- "RNA"
parotid <- SplitObject(parotid, split.by = "Treatment")


SMG$PG_ctrl <- parotid$CTRL
SMG$PG_ir <- parotid$IR



# NORMALIZE
SMG <- lapply(X = SMG, FUN = function(y){
  y <- SCTransform(y, vars.to.regress = "percent.mt", verbose = T)
})

features <- SelectIntegrationFeatures(SMG, nfeatures = 3000)
SMG <- PrepSCTIntegration(SMG, anchor.features = features)
anchors <- FindIntegrationAnchors(SMG, anchor.features = features, normalization.method = "SCT")
SMG <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = 70)
SMG <- RunPCA(SMG, npcs = 30)
SMG <- RunUMAP(SMG, reduction = "pca", dims = 1:30)
SMG <- FindNeighbors(SMG, dims = 1:30, force.recalc = T)
SMG <- FindClusters(SMG, resolution = 0.6)


DimPlot(SMG, group.by = "Gland")
DimPlot(SMG, group.by = "CellType")

epithelial_cells <- c("Acinar", "Ascl3+ duct", "Basal duct", "Bpifa2+", "Bpifa2+ Proacinar",
                      "Etv1+", "GCT", "Intercalated duct", "Krt19+ duct", "Myoepithelial", 
                      "Smgc+", "Smgc+ Proacinar", "Striated duct")


epi_subset <- SMG[, SMG$CellType %in% epithelial_cells]

library(RColorBrewer)
mypalette <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set3"))
mypalette[10] <-"goldenrod2"

DimPlot(epi_subset, group.by = "CellType", cols =mypalette, pt.size = 1)

epi_subset <- SetIdent(epi_subset, value = "CellType")
DimPlot(epi_subset, cells.highlight = WhichCells(epi_subset, idents = "Smgc+ Proacinar"))


DimPlot(epi_subset, group.by = "integrated_snn_res.0.6", pt.size = 1, label = T, label.size = 6) +NoLegend()
VlnPlot(epi_subset, features = "Etv1", pt.size = 0, group.by = "integrated_snn_res.0.6")


epi_subset<- SetIdent(epi_subset, value = "integrated_snn_res.0.6")
cluster10_sub <- epi_subset[,epi_subset$integrated_snn_res.0.6 %in% c(10, 13)]

DimPlot(cluster10_sub, group.by = "CellType", cols = mypalette)

counts <- data.frame(table(cluster10_sub[,cluster10_sub$integrated_snn_res.0.6==10]$CellType))
counts2 <- data.frame(table(cluster10_sub[,cluster10_sub$integrated_snn_res.0.6==13]$CellType))
counts$cluster <- "cluster10"
counts2$cluster <- "cluster13"
counts <- rbind(counts, counts2)
names(counts)[1:2] <- c("CellType", "Number")

# library
library(ggplot2)

# Grouped
ggplot(counts, aes(fill=CellType, y=Number, x=cluster)) + 
  geom_bar(position="fill", stat="identity") + theme_classic() + 
  scale_fill_manual(values=mypalette)


cluster10 <- FindMarkers(epi_subset, ident.1 = 10, only.pos = T, group.by = "integrated_snn_res.0.6")
cluster13 <- FindMarkers(epi_subset, ident.1 = 13, only.pos = T, group.by = "integrated_snn_res.0.6")

c10.top <- cluster10 %>% top_n(20, avg_log2FC)
c13.top <- cluster13 %>% top_n(20, avg_log2FC)
genes <- c(rownames(c10.top), rownames(c13.top))

epi_subset<- SetIdent(epi_subset, value = "integrated_snn_res.0.6")
matrix <- AverageExpression(epi_subset, features = unique(genes), assays = "SCT")
matrix <- matrix$SCT

pheatmap::pheatmap(matrix, scale = "row", fontsize = 12, cluster_rows = T, cluster_cols = T,
                   border_color = "black")



# TCELL EXHAUSTION EXPLORATION #
immcells <- par.combined[,par.combined$CellType %in% c("T-cells CD4+", "T-cells CD8+",
                                                       "T-cells CD4+CD8+", "T-cells Cxcr6+",
                                                       "T-cells FoxP3+")]

VlnPlot(immcells, features = c('Cd44', 'Sell', 'Il7r', 'Cxcr3'), ncol = 4, split.by = "Treatment")
VlnPlot(immcells, features = c('Pdcd1', 'Lag3', 'Il7', 'Cxcr3'), ncol = 4, split.by = "Treatment")
