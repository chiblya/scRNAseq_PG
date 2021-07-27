library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(Hmisc)
library(circlize)
library(Matrix)
library(viridis)

## assign custom colors to cell types for visualization
celltype.colors <- c(
  "Acinar 1" ="blue2",
  "Acinar 2" = "dodgerblue2",
  "Intercalated duct" = "cyan2",
  "Striated duct" = "aquamarine3",
  "Myoepithelial" = "firebrick2",
  "Stromal" = "purple",
  "Endothelial" ="pink",
  "B-cells" = "gray",
  "T-cells CD4+" = "moccasin",
  "T-cells CD4+CD8+" = "navajowhite2",
  "T-cells CD8+" = "navajowhite3",
  "T-cells FoxP3+" = "navajowhite4",
  "T-cells Cxcr6+" = "orange",
  "M0" = "olivedrab1",
  "Dendritic cells" = "olivedrab4",
  "NK cells" = "black")
celltype.colors2 <- c(
  "Acinar 1" ="blue2",
  "Acinar 2" = "dodgerblue2",
  "Intercalated duct" = "cyan2",
  "Striated duct" = "aquamarine3",
  "Myoepithelial" = "firebrick2",
  "Stromal" = "purple",
  "Endothelial" ="pink",
  "Immune" = "gray")

celltype.levels <- c("Acinar 1", "Acinar 2", "Intercalated duct", "Striated duct", "Myoepithelial", "Stromal", "Endothelial", "B-cells", "T-cells CD4+", "T-cells CD4+CD8+", "T-cells CD8+", "T-cells FoxP3+", "T-cells Cxcr6+", "M0", "Dendritic cells", "NK cells")
celltype.levels2 <- c("Acinar 1", "Acinar 2", "Intercalated duct", "Striated duct", "Myoepithelial", "Stromal", "Endothelial", "Immune")

### Load data
par.UT.data <- Read10X(data.dir = "~/Alex Manuscripts/MEC manuscript/CellRanger Files/NRTN (2nd run)/NonIR-PAR/filtered_feature_bc_matrix/")
par.IR.data <- Read10X(data.dir = "~/Alex Manuscripts/MEC manuscript/CellRanger Files/NRTN (2nd run)/IR-PAR/filtered_feature_bc_matrix/")
ut.par <- CreateSeuratObject(counts = par.UT.data, min.cells = 3, min.features = 100, project = "UTvsIR")
ir.par <- CreateSeuratObject(counts = par.IR.data, min.cells = 3,min.features = 100, project = "UTvsIR")

ut.par[["Treatment"]] <- "CTRL"
ir.par[["Treatment"]] <- "IR"
ut.par[["Gland"]] <- "Parotid"
ir.par[["Gland"]] <- "Parotid"

### Normalize and scale using standard worflow
### we create a list to reduce code and simplify process
### list ensures both datasets are normalized and scaled using same method

my.list <- list(ut.par,ir.par)
for (i in 1:length(my.list)) {
  my.list[[i]][["percent.mt"]] <- PercentageFeatureSet(my.list[[i]], pattern = "^mt-")
  my.list[[i]] <- subset(my.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 & nCount_RNA < 15000)
  my.list[[i]] <- NormalizeData(my.list[[i]], verbose = FALSE)
  my.list[[i]] <- FindVariableFeatures(my.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  all.genes<- rownames(my.list[[i]])
  my.list[[i]] <- ScaleData(my.list[[i]], features = all.genes, vars.to.regress = c("percent.mt"))
  my.list[[i]] <- RunPCA(my.list[[i]], npcs = 30, verbose = FALSE)
}

ElbowPlot(my.list[[1]])
ElbowPlot(my.list[[2]])

## we use clustree to determine optimal resolution for clustering
## optimal resolution will be determined for control sample and
## then both will be clustered at the same resolution for consistency

my.list[[1]] <- FindNeighbors(my.list[[1]], dims = 1:12)
for (resolution in c(0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.4))
{my.list[[1]] <- FindClusters(my.list[[1]], resolution = resolution)
}
library(clustree)
# Perform non-linear dimensional reduction
pdf('clustree_seurat.pdf', width = 30, height = 20)
clustree(my.list[[1]])
dev.off()

### Optimal resolution is at around 0.5 but a preliminary analysis suggests that
### there are known supopulations within clusters that aren't properly separated
### we'll take a different approach only for annotating clusters:
### first, we increase the resolution but label clusters with a broad annotation
### based on lineage (epithelial, immune, etc). Then each individual lineage
### is subclustered only to visualize different cell populations and assign proper labels

for (i in 1:length(my.list)){
  my.list[[i]] <- FindNeighbors(my.list[[i]], dims = 1:12)
  my.list[[i]] <- FindClusters(my.list[[i]], resolution = 0.8)
  my.list[[i]] <- RunUMAP(my.list[[i]], dims = 1:12)
}

### extract individual control and IR files from list
par.control <- my.list[[1]]
par.ir <- my.list[[2]]

###remove unnecessary files to save space
remove(ut.par, ir.par, my.list)

### Unsupervised cluster tSNE plots
DimPlot(par.control, pt.size =1, label = T, label.size = 5) + NoLegend()
DimPlot(par.ir, pt.size =1, label = T, label.size = 5) + NoLegend()

### Identify cluster markers (unsupervised)
control.markers <- FindAllMarkers(par.control, only.pos = T, verbose = T)
control.markers <- control.markers[control.markers$p_val_adj <0.05, ]
ir.markers <- FindAllMarkers(par.ir, only.pos = T, verbose = T)
ir.markers <- ir.markers[ir.markers$p_val_adj<0.05, ]

control.markers.top <- control.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
ir.markers.top <- ir.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

## Rename clusters based on a broad lineage identity
## this identity will be used to subcluster epithelial and immune cells (major populations)
## to assign proper labels to individual cell types based on current knowledge of SG cell types
par.control <- SetIdent(par.control, value = "seurat_clusters")
FeaturePlot(par.control, "Col3a1") # just to confirm presence of stromal cells (cluster 14)
par.control <- RenameIdents(par.control, '0' = "Immune", '1'= "Immune", '2' = "Epithelial", '3' = "Immune", '4'= "Immune", '5' = "Immune", '6'=  "Immune", '7'= "Epithelial", '8' = "Epithelial", '9' = "Immune", '10' = "Epithelial", '11' = "Endothelial", '12'=  "Immune", '13' = "Immune", '14' = "Stromal", '15'="Immune")
par.control[["broad.annotation"]] <- Idents(par.control)

### separate immune populations to recluster and annotate
par.immune <- subset(par.control, idents = "Immune")
par.immune <- NormalizeData(par.immune)
par.immune <- FindVariableFeatures(par.immune, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
all.genes <- rownames(par.immune)
par.immune <- ScaleData(par.immune, features = all.genes)
par.immune <- RunPCA(par.immune, npcs = 30, verbose = FALSE)
ElbowPlot(par.immune)

par.immune <- FindNeighbors(par.immune, dims = 1:10)
par.immune <- FindClusters(par.immune, resolution = 0.8)
par.immune <- RunUMAP(par.immune, reduction = "pca", dims = 1:10)

DimPlot(par.immune, group.by = "seurat_clusters", repel = F,pt.size = 1,label = T,label.size = 6) + theme(axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 16, color = "black"))

immune.par.markers <- FindAllMarkers(par.immune, logfc.threshold = 0.3, min.pct = 0.25, only.pos = T)
immune.par.markers <- immune.par.markers[immune.par.markers$p_val_adj<0.05, ]
immune.par.markers.top <- immune.par.markers%>% group_by(cluster) %>% top_n(10, avg_logFC)

FeaturePlot(par.immune, "Flt3") #marker for B-cells
par.immune <- SetIdent(par.immune, value = "seurat_clusters")
par.immune <- RenameIdents(par.immune, '0' = "T-cells CD4+", '1'= "B-cells", 
                           '2' = "B-cells", '3' = "T-cells FoxP3+", 
                           '4'= "T-cells CD4+CD8+", '5' = "T-cells CD8+",
                           '6'=  "M0", '7'= "B-cells", 
                           '8' = "T-cells Cxcr6+", '9'="M0", '10'="NK cells", 
                           '11'="Dendritic cells")
# DCells <- CellSelector(FeaturePlot(par.immune, features = "Flt3"))
# Idents(par.immune, cells=DCells) <- "Dendritic cells"
par.immune[["CellType.fixed"]] <-Idents(par.immune)
DimPlot(par.immune, label = T, repel = T) +NoLegend()

# transfer labels to main seurat object with all populations
select.cells <- list() 
for (i in 1:9) {
  select.cells[[i]] <- WhichCells(object = par.immune, idents = levels(Idents(par.immune))[[i]])
  par.control <- SetIdent(par.control, cells = select.cells[[i]],value = levels(Idents(par.immune))[[i]])
}
DimPlot(par.control)

#### repeat process to label epithelial clusters
par.epithelium <- subset(par.control, idents = c("Epithelial"))
par.epithelium <- NormalizeData(par.epithelium)
par.epithelium <- FindVariableFeatures(par.epithelium, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
all.genes <- rownames(par.epithelium)
par.epithelium <- ScaleData(par.epithelium, features = all.genes)
par.epithelium <- RunPCA(par.epithelium, npcs = 30, verbose = FALSE)
ElbowPlot(par.epithelium)
par.epithelium <- FindNeighbors(par.epithelium, reduction = "pca", dims = 1:8, force.recalc = T)
par.epithelium <- FindClusters(par.epithelium, resolution = 0.8)
par.epithelium <- RunUMAP(par.epithelium, reduction = "pca", dims = 1:8)

DimPlot(par.epithelium, group.by = "seurat_clusters",pt.size = 1,label = T,label.size = 6) + theme(axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 16, color = "black")) + NoLegend()

epi.par.markers <- FindAllMarkers(par.epithelium, logfc.threshold = 0.3, min.pct = 0.25, only.pos = T)
epi.par.markers <- epi.par.markers[epi.par.markers$p_val_adj<0.05, ]
epipar.markers.top <- epi.par.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

FeaturePlot(par.epithelium, features = "Etv1", min.cutoff = "q10", max.cutoff = "q90")
par.epithelium <- SetIdent(par.epithelium, value = "seurat_clusters")
par.epithelium <- RenameIdents(par.epithelium, '0' = "Acinar 1", '1'= "Acinar 1",
                               '2' = "Acinar 1", '3' = "Striated duct", '4'= "Immune", 
                               '5' = "Intercalated duct", '6'=  "Acinar 2", '7'="Immune")
## clusters 4 and 7 were comprised of cells with both immune and epithelial markers
## and are likely doublets. We maintain the label "immune" to visualize them later 

par.epithelium[["CellType.fixed"]] <-Idents(par.epithelium)
DimPlot(par.epithelium)

select.cells <- list()
for (i in 1:5) {
  select.cells[[i]] <- WhichCells(object = par.epithelium, idents = levels(Idents(par.epithelium))[[i]])
  par.control <- SetIdent(par.control, cells = select.cells[[i]],value = levels(Idents(par.epithelium))[[i]])
}

## we can see with these plots that the cells labeled "immune" which express epithelial and immune markers
## and are distributed throughout the UMAP space and are not delimited to a defined cluster.
## These are likely doublets and will be discarded

DimPlot(par.control, cells.highlight = WhichCells(object = par.epithelium,idents = "Immune"))
VlnPlot(par.epithelium, features = "percent.mt")

DimPlot(par.control)
par.control[["CellType.fixed"]] <- Idents(par.control)

par.control.filtered <- subset(par.control, idents = "Immune", invert=T)
DimPlot(par.control.filtered)


par.control.filtered <- NormalizeData(par.control.filtered)
par.control.filtered <- FindVariableFeatures(par.control.filtered, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
all.genes <- rownames(par.control.filtered)
par.control.filtered <- ScaleData(par.control.filtered, features = all.genes, vars.to.regress = "percent.mt")
par.control.filtered <- RunPCA(par.control.filtered, npcs = 30, verbose = FALSE)
ElbowPlot(par.control.filtered)
par.control.filtered <- FindNeighbors(par.control.filtered, reduction = "pca", dims = 1:12)
par.control.filtered <- FindClusters(par.control.filtered, resolution = 0.8)
par.control.filtered <- RunTSNE(par.control.filtered, reduction = "pca", dims = 1:12)
DimPlot(par.control.filtered, group.by = "CellType.fixed", label = T)

par.control.filtered <- RenameIdents(par.control.filtered, '0' = "B-cells", 
                                     '1'= "T-cells CD4+CD8+", '2' = "Acinar 1",
                                     '3' = "B-cells", '4'= "T-cells CD4+", 
                                     '5' = "T-cells CD8+", '6'=  "T-cells FoxP3+",
                                     '7'= "M0", '8' = "Acinar 2", '9' = "Striated duct",
                                     '10' = "Intercalated duct", '11' = "Endothelial", 
                                     '12'=  "T-cells Cxcr6+", '13' = "Dendritic cells", 
                                     '14' = "Stromal", '15'="Dendritic cells")

### NK Cells and Myoepithelial were manually annotated
MECs <- CellSelector(FeaturePlot(par.control.filtered,features = "Acta2" ))
NKcells <- CellSelector(FeaturePlot(par.control.filtered, "Gzma"))
Idents(par.control.filtered, cells=MECs) <- "Myoepithelial"
Idents(par.control.filtered, cells=NKcells) <- "NK cells"
Idents(par.control.filtered) <- factor(Idents(par.control.filtered), levels = celltype.levels)

### Dimplot for figure
DimPlot(par.control.filtered, cols = celltype.colors, pt.size = 1)
DimPlot(par.control.filtered, group.by = "seurat_clusters", pt.size = 1, label = T, label.size = 5)


## IDENTIFY CELL MARKERS
par.control.filtered$CellType.fixed <- Idents(par.control.filtered)
control.markers <- FindAllMarkers(par.control.filtered, only.pos = T)
control.markers <- control.markers[control.markers$p_val_adj <0.05, ]
control.markers.top <- control.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
write.csv(control.markers, file = "Control PG Cell Type Markers.csv")


par.control.filtered<- SetIdent(par.control.filtered,value = "seurat_clusters")
unsup.markers <- FindAllMarkers(par.control.filtered, only.pos = T)
unsup.markers <- unsup.markers[unsup.markers$p_val_adj <0.05, ]
unsup.markers.top <- unsup.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
write.csv(unsup.markers, file = "Control PG Unsupervised Markers.csv")

## DOT PLOT FIGURE 1
DotPlot(par.control.filtered, features = rev(unique(control.markers.top$gene[1:80])),col.min = 0, cols = "Spectral", dot.min = 0.05, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
DotPlot(par.control.filtered, features = rev(unique(unsup.markers.top$gene[1:80])),col.min = 0, cols = "Spectral", dot.min = 0.05) + theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

### Table with count number of cells per cluster
cells.per.cluster.unsup <- as.data.frame(table(par.control.filtered$seurat_clusters))
cells.per.cluster.celltype <-as.data.frame(table(par.control.filtered$CellType.fixed))

saveRDS(par.control.filtered, file = "PG control annotated (Fig1).rds")

#### FIGURE 2 -- ACINAR 2 POPULATION #######
par.control.filtered<- SetIdent(par.control.filtered,value = "seurat_clusters")

DimPlot(subset(par.control.filtered, idents = c(2,8,9,10)), pt.size = 1, label = F, label.size = 6) + NoLegend()
DimPlot(par.control.filtered, group.by = "Treatment",cols = "gray", pt.size = 1)
genes.to.heatmap <- unsup.markers[unsup.markers$cluster %in% c("2", "8", "9", "10"),]$gene 
DoHeatmap(subset(par.control.filtered,idents = c(2,8,9,10)),features = genes.to.heatmap, disp.min = -1 ) + theme(axis.text = element_blank()) + scale_fill_viridis(option = "magma")


par.control.filtered<- SetIdent(par.control.filtered,value = "CellType.fixed")
acinar1markers <- control.markers[control.markers$cluster %in% "Acinar 1",]
acinar1markers <- acinar1markers[order(acinar1markers$avg_logFC, decreasing = T), ]
acinar1markers <- acinar1markers[acinar1markers$p_val_adj<0.05, ]
  
acinar2markers <- control.markers[control.markers$cluster %in% "Acinar 2",]
acinar2markers <- acinar2markers[order(acinar2markers$avg_logFC, decreasing = T), ]
acinar2markers <- acinar2markers[acinar2markers$p_val_adj<0.05, ]

common.acinar.markers <- acinar1markers[acinar1markers$gene %in% acinar2markers$gene, ]
acinar1.unique.markers <- acinar1markers[!(acinar1markers$gene %in% acinar2markers$gene), ]
acinar2.unique.markers <- acinar2markers[!(acinar2markers$gene %in% acinar1markers$gene), ]


DotPlot(par.control.filtered, features = rev(unique(acinar2markers$gene[1:50])),col.min = 0, cols = "Spectral", dot.min = 0.05) + theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

acinar1.vs.acinar2 <- FindMarkers(par.control.filtered, ident.1 = "Acinar 1", ident.2 = "Acinar 2", only.pos = F)
acinar1.vs.acinar2 <- acinar1.vs.acinar2[acinar1.vs.acinar2$p_val_adj<0.05,]
acinar1.vs.acinar2$gene <- rownames(acinar1.vs.acinar2)
acinar1.vs.acinar2 <- acinar1.vs.acinar2[order(acinar1.vs.acinar2$avg_logFC), ]

VlnPlot(subset(par.control.filtered, idents = c("Acinar 1", "Acinar 2")), features = tail(acinar1.vs.acinar2$gene, 10), pt.size = 0, ncol  = 10)
VlnPlot(subset(par.control.filtered, idents = c("Acinar 1", "Acinar 2")), features = c("Amy1", "Bpifa2"), pt.size = 0, ncol  = 2)
VlnPlot(subset(par.control.filtered, idents = c("Acinar 1", "Acinar 2")), features = c("Pax9", "Atf3"), pt.size = 0, ncol  = 2)


########### LIGAND RECEPTOR ANALYSIS ##############
##############################################################
########## Ligand-Receptor analysis ##########################

### Load Ligand-Receptor Pair info
Ligand.Receptor.Pairs <- read.delim("~/Alex Manuscripts/MEC manuscript/Figure 2 - MEC Functions/Ligand-Receptor Pairs.txt")
Ligand.Receptor.Pairs$Pair <- tolower(Ligand.Receptor.Pairs$Pair)
Ligand.Receptor.Pairs$Pair <- capitalize(Ligand.Receptor.Pairs$Pair)
Ligand.Receptor.Pairs$Ligand <- tolower(Ligand.Receptor.Pairs$Ligand)
Ligand.Receptor.Pairs$Ligand <- capitalize(Ligand.Receptor.Pairs$Ligand)
Ligand.Receptor.Pairs$Receptor <- tolower(Ligand.Receptor.Pairs$Receptor)
Ligand.Receptor.Pairs$Receptor <- capitalize(Ligand.Receptor.Pairs$Receptor)

Ligand.Receptor.Pairs <- data.frame(lapply(Ligand.Receptor.Pairs, function(x) {
  gsub("Ntf4", "Ntf5", x) 
  }))


##### Create simple cell type labels for better visualization #####

par.control.filtered <- SetIdent(par.control.filtered, value = "CellType.fixed")
par.control.filtered <- RenameIdents(par.control.filtered, "B-cells" = "Immune",
                                     "Dendritic cells"="Immune",
                                     "M0"="Immune",
                                     "NK cells"="Immune",
                                     "T-cells CD4+"="Immune",
                                     "T-cells CD4+CD8+"="Immune",
                                     "T-cells CD8+"="Immune",
                                     "T-cells Cxcr6+"="Immune",
                                     "T-cells FoxP3+"="Immune")

identities <- as.vector(unique(Idents(par.control.filtered)))
identities <- factor(identities, levels = celltype.levels2)

Idents(par.control.filtered) <- factor(Idents(par.control.filtered), levels = celltype.levels2 )
par.control.filtered$CellType.simple <- Idents(par.control.filtered) 
DimPlot(par.control.filtered)

simple.control.markers <- FindAllMarkers(par.control.filtered, only.pos = T)
simple.control.markers <- simple.control.markers[simple.control.markers$p_val_adj<0.05, ]

##### count highly expressed ligands and receptors only (defining genes) ###
ligand_list <- vector(mode = "list", length = 8)
names(ligand_list) <- names(celltype.colors2)

for (i in 1:8) {
  temptable <- simple.control.markers[simple.control.markers$cluster %in% names(ligand_list)[i],]
  ligand_list[[i]] <- temptable[temptable$gene %in% Ligand.Receptor.Pairs$Ligand, ]
}

receptor_list <- vector(mode = "list", length = 8)
names(receptor_list) <- names(celltype.colors2)

for (i in 1:8) {
  temptable <- simple.control.markers[simple.control.markers$cluster %in% names(receptor_list)[i],]
  receptor_list[[i]] <- temptable[temptable$gene %in% Ligand.Receptor.Pairs$Receptor, ]
}

receptor.counter = c(0,0,0,0,0,0,0,0)
ligand.counter = c(0,0,0,0,0,0,0,0)
names(receptor.counter) <- names(ligand_list)
names(ligand.counter) <- names(ligand_list)
for (i in 1:8) {
  receptor.counter[i]<- dim(receptor_list[[i]])[1]
  ligand.counter[i]<- dim(ligand_list[[i]])[1]
}
receptor.counter
ligand.counter

par(mar=c(8,4,4,1))
barplot(receptor.counter,
        ylab="Number of enriched receptors",
        border="black", col = "blue", las=2, axis.lty = 1)

barplot(ligand.counter,
        ylab="Number of enriched ligands",
        border="black", col = "blue", las=2, axis.lty = 1
)

#### Make table with Strong putative ligand-receptor pairs
#### based on cluster-defining genes

library(reshape2)
allfoldchanges <- dcast(simple.control.markers,formula = gene~cluster,fun.aggregate = sum,value.var = "avg_logFC")
AllFC.receptors <- merge(allfoldchanges, Ligand.Receptor.Pairs, by.x = "gene", by.y = "Receptor")
AllFC.Ligands <- merge(allfoldchanges, Ligand.Receptor.Pairs, by.x = "gene", by.y = "Ligand")

Strong.potential.Pairs <- merge(AllFC.receptors, allfoldchanges, by.x = "Ligand", by.y = "gene", no.dups = F)
names(Strong.potential.Pairs)[2] <- "Receptor"
write.table(Strong.potential.Pairs, file = "Strong potential Ligand-Receptor Pairs (based on cluster-defining genes only).txt", sep = "\t", row.names = T)


###### ACINAR LIGANDS AND RECEPTORS
DotPlot(par.control.filtered, features = rev(ligand_list$`Acinar 1`$gene), dot.min = 0.05,dot.scale = 5) +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(par.control.filtered, features = rev(ligand_list$`Acinar 2`$gene), dot.min = 0.05,dot.scale = 5) +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(par.control.filtered, features = rev(receptor_list$`Acinar 1`$gene),scale.max = 100, scale.min = 5) +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(par.control.filtered, features = rev(receptor_list$`Acinar 2`$gene), dot.min = 0.05,dot.scale = 5) +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


##### create matrix with number of potential pairs
counter = c(0,0,0,0,0,0,0,0)  #create counter vector with 8 slots for 8 cell types
number.of.pairs <- matrix(ncol = 8, nrow = 8, dimnames = list(colnames(Strong.potential.Pairs[3:10]), colnames(Strong.potential.Pairs[3:10]) )) # create 10X10 matrix to save number of pairs between cell types
for (j in 1:8){
  counter = c(0,0,0,0,0,0,0,0)
  for (i in 1:NROW(Strong.potential.Pairs)){
    if(Strong.potential.Pairs[i,j+2]>0){
      for(k in 1:8){
        if(Strong.potential.Pairs[i,k+13]>0){
          counter[k] <- counter[k] + 1
        }
      }
    }
  }
  number.of.pairs[j,] <- counter[1:8]
}

#### The transposed matrix of 'Number of Pairs' is organized to plot the recipient (receptor) cells in columns and provider (ligands) cells in rows
#### The resulting plot is ligands from X population towards receptors in all other cells
number.of.pairs <- t(number.of.pairs)

#### chord diagrams
colnames(number.of.pairs) <-  sub(".x", "", colnames(number.of.pairs))
rownames(number.of.pairs) <- sub(".x", "", rownames(number.of.pairs))
write.csv(number.of.pairs, file = "PG Control-Number of ligand-receptor potential pairs with DEGs.csv")

border_mat = matrix("black", nrow = 2, ncol = ncol(number.of.pairs))
rownames(border_mat) = c("Acinar 1",  "Acinar 2")
colnames(border_mat) = colnames(number.of.pairs)
grid.col = celltype.colors2

col_mat = rand_color(length(number.of.pairs), transparency = 0.5)
dim(col_mat) = dim(number.of.pairs)  # to make sure it is a matrix

col_mat[c(1:8),c(1:8)]<- '#00000010'
col_mat[1,] <- celltype.colors2
col_mat[2,] <- "dodgerblue"

dev.off()
circos.par(gap.degree=5, gap.after=5)
chordDiagram(number.of.pairs[1:2,], grid.col = grid.col, link.lwd = 1, link.lty = 1, link.border = border_mat, 
             symmetric = F, directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.width = 0.1, 
             link.arr.length = 0.1, link.arr.type = "big.arrow",
             annotationTrack = "grid", link.largest.ontop = T, link.arr.lty = 1,grid.border = 1,
             row.col = c("#0000FF", "#3399FF", "#00FFFF33","#00999933",
                         "#FF000033", "#7F00FF33", "#FF66FF33", "#E0E0E033"))
circos.clear()

circos.par(gap.degree=5, gap.after=5)
chordDiagram(number.of.pairs, grid.col = grid.col, link.lwd = 1, link.lty = 1, link.border = border_mat, 
             symmetric = F, directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.width = 0.1, 
             link.arr.length = 0.1, link.arr.type = "big.arrow",
             annotationTrack = "grid", link.largest.ontop = T, link.arr.lty = 1,grid.border = 1,
             row.col = col_mat)
circos.clear()

### Plot ligands to acinar populations ###
border_mat = matrix("black", nrow = nrow(number.of.pairs), ncol = 2)
rownames(border_mat) = rownames(number.of.pairs)
colnames(border_mat) =  c("Acinar 1",  "Acinar 2")
grid.col = celltype.colors2

col_mat = rand_color(length(number.of.pairs), transparency = 0.5)
dim(col_mat) = dim(number.of.pairs)  # to make sure it is a matrix

col_mat[c(1:8),c(1:8)]<- '#00000010'
col_mat[,1] <- celltype.colors2
col_mat[,2] <- celltype.colors2

dev.off()
circos.par(gap.degree=5, gap.after=5)
chordDiagram(number.of.pairs[,1:2], grid.col = grid.col, link.lwd = 1, link.lty = 1, link.border = border_mat, 
             symmetric = F, directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.width = 0.1, 
             link.arr.length = 0.1, link.arr.type = "big.arrow",
             annotationTrack = "grid", link.largest.ontop = T, link.arr.lty = 1,grid.border = 1,
             col = col_mat  )
circos.clear()

dev.off()
circos.par(gap.degree=5, gap.after=5)
chordDiagram(number.of.pairs, grid.col = grid.col, link.lwd = 1, link.lty = 1, link.border = border_mat, 
             symmetric = F, directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.width = 0.1, 
             link.arr.length = 0.1, link.arr.type = "big.arrow",
             annotationTrack = "grid", link.largest.ontop = T, link.arr.lty = 1,grid.border = 1,
             col = col_mat  )
circos.clear()

#####
#### ALL EXPRESSED LIGANDS AND RECEPTORS BASED ON EXPRESSION THRESHOLD OF O.1
#### GENES ARE NOT EXLUSIVE TO A CELL CLUSTER

AvgExpression.allclusters <- AverageExpression(par.control.filtered, assays = "RNA")
head(AvgExpression.allclusters[["RNA"]][, 1:5])
allaverages <- as.data.frame(AvgExpression.allclusters[["RNA"]])
write.table(x = allaverages, file = "PG Average gene expression per cluster (simpler annotations).txt", sep = "\t", row.names = T)
indices = abs(allaverages)<0.1
allaverages[indices] = NA
allaverages <- allaverages[!(rowSums(is.na(allaverages))==NCOL(allaverages)),] ##remove rows that contain NA all across
allaverages$gene <- rownames(allaverages)

### Create list of ligands and receptors with expression values per cell type
AllAvgs.UTvsIR.Receptors <- merge(allaverages, Ligand.Receptor.Pairs, by.x = "gene", by.y = "Receptor")
AllAvgs.UTvsIR.Ligands <- merge(allaverages, Ligand.Receptor.Pairs, by.x = "gene", by.y = "Ligand")

### Count number of ligands per cell
count.ligands <- AllAvgs.UTvsIR.Ligands[!duplicated(AllAvgs.UTvsIR.Ligands$gene), ]
counter = c(0,0,0,0,0,0,0,0)
names(counter) <- colnames(count.ligands[,2:9])
for (i in 2:9) {
  counter[i-1]<- sum(!(is.na(count.ligands[,i])))
}
counter

par(mar=c(8,4,4,2))
barplot(counter,
        ylab="Number of expressed ligands",
        border="black", col = "blue", las=2, axis.lty = 1
)

### Count number of receptors per cell
count.receptors <- AllAvgs.UTvsIR.Receptors[!duplicated(AllAvgs.UTvsIR.Receptors$gene), ]
counter = c(0,0,0,0,0,0,0,0)
names(counter) <- colnames(count.receptors[,2:9])
for (i in 2:9) {
  counter[i-1]<- sum(!(is.na(count.receptors[,i])))
}
counter

par(mar=c(8,4,4,2))
barplot(counter,
        ylab="Number of expressed receptors",
        border="black", col = "blue", las=2, axis.lty = 1
)

##### Add ligand expression info to receptor table to identify all potential ligand-receptor pairs
AvgExpression.potential.Pairs <- merge(AllAvgs.UTvsIR.Receptors, allaverages, by.x = "Ligand", by.y = "gene", no.dups = F)
names(AvgExpression.potential.Pairs)[2] <- "Receptor"
write.table(AvgExpression.potential.Pairs, file = "All potential Ligand-Receptor Pairs.txt", sep = "\t", row.names = T)



##### create matrix with number of potential pairs
counter = c(0,0,0,0,0,0,0,0)  #create counter vector with 8 slots for 8 cell types
number.of.pairs <- matrix(ncol = 8, nrow = 8, dimnames = list(colnames(Strong.potential.Pairs[3:10]), colnames(Strong.potential.Pairs[3:10]) )) # create 10X10 matrix to save number of pairs between cell types
for (j in 1:8){
  counter = c(0,0,0,0,0,0,0,0)
  for (i in 1:NROW(Strong.potential.Pairs)){
    if(Strong.potential.Pairs[i,j+2]>0){
      for(k in 1:8){
        if(Strong.potential.Pairs[i,k+13]>0){
          counter[k] <- counter[k] + 1
        }
      }
    }
  }
  number.of.pairs[j,] <- counter[1:8]
}

#### The transposed matrix of 'Number of Pairs' is organized to plot the recipient (receptor) cells in columns and provider (ligands) cells in rows
#### The resulting plot is ligands from X population towards receptors in all other cells
number.of.pairs <- t(number.of.pairs)

#### chord diagrams
colnames(number.of.pairs) <-  sub(".x", "", colnames(number.of.pairs))
rownames(number.of.pairs) <- sub(".x", "", rownames(number.of.pairs))
border_mat = matrix("black", nrow = 2, ncol = ncol(number.of.pairs))
rownames(border_mat) = c("Acinar 1",  "Acinar 2")
colnames(border_mat) = colnames(number.of.pairs)
grid.col = celltype.colors2

dev.off()
circos.par(gap.degree=5, gap.after=5)
chordDiagram(number.of.pairs, grid.col = grid.col, link.lwd = 1, link.lty = 1, link.border = border_mat, 
             symmetric = F, directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.width = 0.1, 
             link.arr.length = 0.1, link.arr.type = "big.arrow",
             annotationTrack = "grid", link.largest.ontop = T, link.arr.lty = 1,grid.border = 1,
             row.col = c("#0000FF", "#3399FF", "#00FFFF33","#00999933",
                         "#FF000033", "#7F00FF33", "#FF66FF33", "#E0E0E033"))
circos.clear()


### Plot ligands to acinar populations ###
border_mat = matrix("black", nrow = nrow(number.of.pairs), ncol = 2)
rownames(border_mat) = rownames(number.of.pairs)
colnames(border_mat) =  c("Acinar 1",  "Acinar 2")
grid.col = celltype.colors2


col_mat = rand_color(length(number.of.pairs), transparency = 0.5)
dim(col_mat) = dim(number.of.pairs)  # to make sure it is a matrix

col_mat[c(1:8),c(1:8)]<- '#00000010'
col_mat[,1] <- celltype.colors2
col_mat[,2] <- celltype.colors2

dev.off()
circos.par(gap.degree=5, gap.after=5)
chordDiagram(number.of.pairs, grid.col = grid.col, link.lwd = 1, link.lty = 1, link.border = border_mat, 
             symmetric = F, directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.width = 0.1, 
             link.arr.length = 0.1, link.arr.type = "big.arrow",
             annotationTrack = "grid", link.largest.ontop = T, link.arr.lty = 1,grid.border = 1,
             col = col_mat  )
circos.clear()

