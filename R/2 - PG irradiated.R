library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(Hmisc)
library(circlize)
library(Matrix)
library(viridis)

#################################
##### IR SAMPLE #########
#This script continues from PG Control analysis
####

par.ir <- SetIdent(par.ir, value = "seurat_clusters")
par.ir <- RenameIdents(par.ir, '0' = "Immune", '1'= "Immune", '2' = "Immune", '3' = "Immune", '4'= "Epithelial", '5' = "Immune", '6'=  "Immune", '7'= "Immune", '8' = "Epithelial", '9' = "Immune", '10' = "Immune", '11' = "Immune", '12'=  "Endothelial", '13' = "Epithelial", '14' = "Immune", '15'="Epithelial", '16'="Immune", '17'="Immune", '18'="Immune")
par.ir[["broad.annotation"]] <- Idents(par.ir)
DimPlot(par.ir)

FeaturePlot(par.ir, features = "Acta2")
Stromal <- CellSelector(FeaturePlot(par.ir, features = "Vim"))
MECs <- CellSelector(FeaturePlot(par.ir, features = "Acta2"))
Idents(par.ir, cells=Stromal) <- "Stromal"
Idents(par.ir, cells=MECs) <- "Myoepithelial"
par.ir[["broad.annotation"]] <- Idents(par.ir)

##### Separate Immune populations to achieve better clustering ######
#### Immune clusters
par.ir <- SetIdent(par.ir, value = "broad.annotation")
par.immune <- subset(par.ir, idents = "Immune")
par.immune <- NormalizeData(par.immune)
par.immune <- FindVariableFeatures(par.immune, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
all.genes <- rownames(par.immune)
par.immune <- ScaleData(par.immune, features = all.genes)
par.immune <- RunPCA(par.immune, npcs = 30, verbose = FALSE)
ElbowPlot(par.immune)

par.immune <- FindNeighbors(par.immune, reduction = "pca", dims = 1:10)
par.immune <- FindClusters(par.immune, resolution = 0.8)
par.immune <- RunTSNE(par.immune, reduction = "pca", dims = 1:10)

FeaturePlot(par.immune, features = "Flt3")
DimPlot(par.immune, group.by = "seurat_clusters",pt.size = 1,label = T,label.size = 6) + theme(axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 16, color = "black"))
DimPlot(par.immune, group.by = "Treatment", pt.size = 0.5,label = F,label.size = 8, cols = c("#00009980", "#FF666633"))+ theme(axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 16, color = "black"))

immune.par.markers <- FindAllMarkers(par.immune, logfc.threshold = 0.3, min.pct = 0.25, only.pos = T)
immune.par.markers <- immune.par.markers[immune.par.markers$p_val_adj<0.05, ]
immune.par.markers.top <- immune.par.markers%>% group_by(cluster) %>% top_n(10, avg_logFC)

par.immune <- RenameIdents(par.immune, '0' = "T-cells CD4+", '1'= "T-cells CD4+CD8+", '2' = "T-cells FoxP3+", '3' = "B-cells", '4'= "T-cells CD4+", '5' = "B-cells",
                           '6'=  "NK cells", '7'= "T-cells CD8+", '8' = "T-cells CD8+", '9'="Epithelial", '10'="B-cells", '12'="T-cells Cxcr6+", '11'="M0",
                           '13' = "NK cells",'14'="Dendritic cells", '15'="Dendritic cells", '16'="M0")

DimPlot(par.ir)

select.cells <- list()
for (i in 1:10) {
  select.cells[[i]] <- WhichCells(object = par.immune, idents = levels(Idents(par.immune))[[i]])
  par.ir <- SetIdent(par.ir, cells = select.cells[[i]],value = levels(Idents(par.immune))[[i]])
}
DimPlot(par.ir)

#### annotate epithelium clusters
par.epithelium <- subset(par.ir, idents = "Epithelial")
par.epithelium <- NormalizeData(par.epithelium)
par.epithelium <- FindVariableFeatures(par.epithelium, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
all.genes <- rownames(par.epithelium)
par.epithelium <- ScaleData(par.epithelium, features = all.genes)
par.epithelium <- RunPCA(par.epithelium, npcs = 30, verbose = FALSE)
ElbowPlot(par.epithelium)
par.epithelium <- FindNeighbors(par.epithelium, reduction = "pca", dims = 1:5, force.recalc = T)
par.epithelium <- FindClusters(par.epithelium, resolution = 0.5)
par.epithelium <- RunTSNE(par.epithelium, reduction = "pca", dims = 1:5)

DimPlot(par.epithelium, group.by = "seurat_clusters",pt.size = 1,label = T,label.size = 6) + theme(axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 16, color = "black")) + NoLegend()
DimPlot(par.epithelium, group.by = "Treatment", pt.size = 0.5,label = F,label.size = 8, cols = c("#00009980", "#FF666633"))+ theme(axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 16, color = "black"))

epi.par.markers <- FindAllMarkers(par.epithelium, logfc.threshold = 0.3, min.pct = 0.25, only.pos = T)
epi.par.markers <- epi.par.markers[epi.par.markers$p_val_adj<0.05, ]
epipar.markers.top <- epi.par.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

par.epithelium <- SetIdent(par.epithelium, value = "seurat_clusters")

par.epithelium <- RenameIdents(par.epithelium, '0' = "Acinar 1", '1'= "Immune", '2' = "Intercalated duct", 
                               '3' = "Immune", '4'="Acinar 2", '5'= "Striated duct", '6' = "Intercalated duct")


### replace identities in par.ir object
select.cells <- list()
for (i in 1:5) {
  select.cells[[i]] <- WhichCells(object = par.epithelium, idents = levels(Idents(par.epithelium))[[i]])
  par.ir <- SetIdent(par.ir, cells = select.cells[[i]],value = levels(Idents(par.epithelium))[[i]])
}

DimPlot(par.ir,cells.highlight = WhichCells(par.ir, idents = "Immune"))
par.ir[["CellType.fixed"]] <- Idents(par.ir)
par.ir.filtered <- subset(par.ir, idents = "Immune", invert=T)
DimPlot(par.ir.filtered, group.by = "CellType.fixed")

par.ir.filtered <- NormalizeData(par.ir.filtered)
par.ir.filtered <- FindVariableFeatures(par.ir.filtered, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
all.genes <- rownames(par.ir.filtered)
par.ir.filtered <- ScaleData(par.ir.filtered, features = all.genes, vars.to.regress = "percent.mt")
par.ir.filtered <- RunPCA(par.ir.filtered, npcs = 30, verbose = FALSE)
ElbowPlot(par.ir.filtered)
par.ir.filtered <- FindNeighbors(par.ir.filtered, reduction = "pca", dims = 1:12)
par.ir.filtered <- FindClusters(par.ir.filtered, resolution = 0.8)
par.ir.filtered <- RunTSNE(par.ir.filtered, reduction = "pca", dims = 1:12)

DimPlot(par.ir.filtered, group.by = "seurat_clusters", label = T)
DimPlot(par.ir.filtered, group.by = "CellType.fixed", label = T)

ir.markers <- FindAllMarkers(par.ir.filtered, logfc.threshold = 0.3, min.pct = 0.25, only.pos = T)
ir.markers <- ir.markers[ir.markers$p_val_adj<0.05, ]
ir.markers.top <- ir.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

FeaturePlot(par.ir.filtered, features = "Fxyd2")
par.ir.filtered <- RenameIdents(par.ir.filtered, '0' = "T-cells CD4+", 
                                '1'= "T-cells CD4+CD8+", '2' = "B-cells",
                                '3' = "T-cells FoxP3+", '4'= "T-cells CD4+CD8+", 
                                '5' = "T-cells CD4+", '6'=  "B-cells",
                                '7'= "Acinar 1", '8' = "T-cells CD8+", '9' = "T-cells CD8+",
                                '10' = "Intercalated duct", '11' = "Acinar 2", '12'="M0",
                                '13' ="Endothelial", '14'= "NK cells",
                                '15'=  "T-cells Cxcr6+", '16' = "Dendritic cells", 
                                '17' = "M0", '18'="Dendritic cells")

Stromal <- CellSelector(FeaturePlot(par.ir.filtered, features = "Vim"))
MECs <- CellSelector(FeaturePlot(par.ir.filtered, features = "Acta2"))
duct <- CellSelector(FeaturePlot(par.ir.filtered, features = "Fxyd2"))

Idents(par.ir.filtered, cells=Stromal) <- "Stromal"
Idents(par.ir.filtered, cells=MECs) <- "Myoepithelial"
Idents(par.ir.filtered, cells=duct) <- "Striated duct"
Idents(par.ir.filtered) <- factor(Idents(par.ir.filtered), levels = celltype.levels)

### Dimplot for figure
par.ir.filtered[["CellType.fixed"]] <- Idents(par.ir.filtered)
DimPlot(par.ir.filtered, cols = celltype.colors, pt.size = 1, group.by = "CellType.fixed")
DimPlot(par.ir.filtered, group.by = "seurat_clusters", pt.size = 1, label = T, label.size = 5)

## IDENTIFY CELL MARKERS
ir.markers <- FindAllMarkers(par.ir.filtered, only.pos = T)
ir.markers <- ir.markers[ir.markers$p_val_adj <0.05, ]
ir.markers.top <- ir.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
write.csv(ir.markers, file = "IR PG Cell Type Markers.csv")

par.ir.filtered<- SetIdent(par.ir.filtered,value = "seurat_clusters")
unsup..ir.markers <- FindAllMarkers(par.ir.filtered, only.pos = T)
unsup..ir.markers <- unsup..ir.markers[unsup..ir.markers$p_val_adj <0.05, ]
unsup..ir.markers.top <- unsup..ir.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
write.csv(unsup..ir.markers, file = "IR PG Unsupervised Markers.csv")

## DOT PLOT FIGURE 4
DotPlot(par.ir.filtered,group.by = "CellType.fixed", features = rev(unique(ir.markers.top$gene[1:80])),col.min = 0, cols = "Spectral", dot.min = 0.05, dot.scale = 5) + theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
DotPlot(par.ir.filtered, features = rev(unique(unsup..ir.markers.top$gene)),col.min = 0, cols = "Spectral", dot.min = 0.05) + theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

### Table with count number of cells per cluster
cells.per.cluster.ir.unsup <- as.data.frame(table(par.ir.filtered$seurat_clusters))
cells.per.cluster.ir.celltype <-as.data.frame(table(par.ir.filtered$CellType.fixed))

saveRDS(par.ir.filtered, file = "PG IR annotated (Fig4).rds")

######## LIGAND RECEPTOR ANALYSIS
##### Create simple cell type labels for better visualization #####

par.ir.filtered <- SetIdent(par.ir.filtered, value = "CellType.fixed")
par.ir.filtered <- RenameIdents(par.ir.filtered, "B-cells" = "Immune",
                                     "Dendritic cells"="Immune",
                                     "M0"="Immune",
                                     "NK cells"="Immune",
                                     "T-cells CD4+"="Immune",
                                     "T-cells CD4+CD8+"="Immune",
                                     "T-cells CD8+"="Immune",
                                     "T-cells Cxcr6+"="Immune",
                                     "T-cells FoxP3+"="Immune")

identities <- as.vector(unique(Idents(par.ir.filtered)))
identities <- factor(identities, levels = celltype.levels2)

Idents(par.ir.filtered) <- factor(Idents(par.ir.filtered), levels = celltype.levels2 )
par.ir.filtered$CellType.simple <- Idents(par.ir.filtered) 
DimPlot(par.ir.filtered)

simple.ir.markers <- FindAllMarkers(par.ir.filtered, only.pos = T)
simple.ir.markers <- simple.ir.markers[simple.ir.markers$p_val_adj<0.05, ]

##### count highly expressed ligands and receptors only (defining genes) ###
ligand_list <- vector(mode = "list", length = 8)
names(ligand_list) <- names(celltype.colors2)

for (i in 1:8) {
  temptable <- simple.ir.markers[simple.ir.markers$cluster %in% names(ligand_list)[i],]
  ligand_list[[i]] <- temptable[temptable$gene %in% Ligand.Receptor.Pairs$Ligand, ]
}

receptor_list <- vector(mode = "list", length = 8)
names(receptor_list) <- names(celltype.colors2)

for (i in 1:8) {
  temptable <- simple.ir.markers[simple.ir.markers$cluster %in% names(receptor_list)[i],]
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
as.data.frame(receptor.counter)
as.data.frame(ligand.counter)

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
allfoldchanges <- dcast(simple.ir.markers,formula = gene~cluster,fun.aggregate = sum,value.var = "avg_logFC")
AllFC.receptors <- merge(allfoldchanges, Ligand.Receptor.Pairs, by.x = "gene", by.y = "Receptor")
AllFC.Ligands <- merge(allfoldchanges, Ligand.Receptor.Pairs, by.x = "gene", by.y = "Ligand")

Strong.potential.Pairs.ir <- merge(AllFC.receptors, allfoldchanges, by.x = "Ligand", by.y = "gene", no.dups = F)
names(Strong.potential.Pairs.ir)[2] <- "Receptor"
write.table(Strong.potential.Pairs.ir, file = "Strong potential Ligand-Receptor Pairs in IR gland (based on cluster-defining genes only).txt", sep = "\t", row.names = T)


###### ACINAR LIGANDS AND RECEPTORS
DotPlot(par.ir.filtered, features = rev(ligand_list$`Acinar 1`$gene), dot.min = 0.05,dot.scale = 5) +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(par.ir.filtered, features = rev(ligand_list$`Acinar 2`$gene), dot.min = 0.05,dot.scale = 5) +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

DotPlot(par.ir.filtered, features = rev(receptor_list$`Acinar 1`$gene),scale.max = 100, scale.min = 5) +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(par.ir.filtered, features = rev(receptor_list$`Acinar 2`$gene), dot.min = 0.05,dot.scale = 5) +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

##### create matrix with number of potential pairs
counter = c(0,0,0,0,0,0,0,0)  #create counter vector with 8 slots for 8 cell types
number.of.pairs <- matrix(ncol = 8, nrow = 8, dimnames = list(colnames(Strong.potential.Pairs.ir[3:10]), colnames(Strong.potential.Pairs.ir[3:10]) )) # create 10X10 matrix to save number of pairs between cell types
for (j in 1:8){
  counter = c(0,0,0,0,0,0,0,0)
  for (i in 1:NROW(Strong.potential.Pairs.ir)){
    if(Strong.potential.Pairs.ir[i,j+2]>0){
      for(k in 1:8){
        if(Strong.potential.Pairs.ir[i,k+13]>0){
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
write.csv(number.of.pairs, file = "PG IR-Number of ligand-receptor potential pairs with DEGs.csv")


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


