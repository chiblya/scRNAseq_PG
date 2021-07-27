library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(Hmisc)
library(circlize)
library(Matrix)
library(viridis)

#################################################################
######### IMMUNE LIGAND RECEPTOR ANALYSIS CTRL VS IR ############
#This script continues from Merged PG CTRL and IR  analysis #####
#################################################################

par.combined.filtered <- SetIdent(par.combined.filtered, value = "CellType.fixed")
DimPlot(par.combined.filtered)
pg.ctrl.ir.markers <- FindAllMarkers(par.combined.filtered, only.pos = T)

library(reshape2)
ctrl.ir.allfoldchanges <- dcast(pg.ctrl.ir.markers,formula = gene~cluster,fun.aggregate = sum,value.var = "avg_logFC")
receptors <- merge(ctrl.ir.allfoldchanges, Ligand.Receptor.Pairs, by.x = "gene", by.y = "Receptor")
ligands <- merge(ctrl.ir.allfoldchanges, Ligand.Receptor.Pairs, by.x = "gene", by.y = "Ligand")

Strong.potential.Pairs.ctrl.ir <- merge(receptors, ctrl.ir.allfoldchanges, by.x = "Ligand", by.y = "gene", no.dups = F)
names(Strong.potential.Pairs.ctrl.ir)[2] <- "Receptor"
colnames(Strong.potential.Pairs.ctrl.ir) <-  sub(".x", "", colnames(Strong.potential.Pairs.ctrl.ir))
names(Strong.potential.Pairs.ctrl.ir)[14] <- "T-cells FoxP3+"
names(Strong.potential.Pairs.ctrl.ir)[15] <- "T-cells Cxcr6+"

write.table(Strong.potential.Pairs.ctrl.ir, file = "Ligand-Receptor Pairs (Ctrl and IR combined).txt", sep = "\t", row.names = T)

##### count highly expressed ligands and receptors only (defining genes) ###
ligand_list.ctrl.ir <- vector(mode = "list", length = 16)
names(ligand_list.ctrl.ir) <- names(celltype.colors)

for (i in 1:16) {
  temptable <- pg.ctrl.ir.markers[pg.ctrl.ir.markers$cluster %in% names(ligand_list.ctrl.ir)[i],]
  ligand_list.ctrl.ir[[i]] <- temptable[temptable$gene %in% Ligand.Receptor.Pairs$Ligand, ]
}

receptor_list.ctrl.ir <- vector(mode = "list", length = 16)
names(receptor_list.ctrl.ir) <- names(celltype.colors)

for (i in 1:16) {
  temptable <- pg.ctrl.ir.markers[pg.ctrl.ir.markers$cluster %in% names(receptor_list.ctrl.ir)[i],]
  receptor_list.ctrl.ir[[i]] <- temptable[temptable$gene %in% Ligand.Receptor.Pairs$Receptor, ]
}

receptor.counter = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
ligand.counter = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
names(receptor.counter) <- names(ligand_list.ctrl.ir)
names(ligand.counter) <- names(ligand_list.ctrl.ir)
for (i in 1:16) {
  receptor.counter[i]<- dim(receptor_list.ctrl.ir[[i]])[1]
  ligand.counter[i]<- dim(ligand_list.ctrl.ir[[i]])[1]
}
as.data.frame(receptor.counter)
as.data.frame(ligand.counter)

par(mar=c(9,4,3,1))
barplot(receptor.counter,
        ylab="Number of enriched receptors",
        border="black", col = "blue", las=2, axis.lty = 1)

barplot(ligand.counter,
        ylab="Number of enriched ligands",
        border="black", col = "blue", las=2, axis.lty = 1
)

###### T cell CD4+ ligands and receptors
DotPlot(par.combined.filtered, features = rev(ligand_list.ctrl.ir$`T-cells CD4+`$gene), dot.min = 0.05,dot.scale = 5) +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(par.combined.filtered, features = rev(receptor_list.ctrl.ir$`T-cells CD4+`$gene),dot.min = 0.05, dot.scale = 5) +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

DotPlot(par.combined.filtered, features = rev(ligand_list.ctrl.ir$`T-cells Cxcr6+`$gene), dot.min = 0.05,dot.scale = 5) +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(par.combined.filtered, features = rev(receptor_list.ctrl.ir$`T-cells Cxcr6+`$gene),dot.min = 0.05,dot.scale = 5) +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

DotPlot(par.combined.filtered, features = rev(ligand_list.ctrl.ir$`Dendritic cells`$gene), dot.min = 0.05,dot.scale = 5) +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(par.combined.filtered, features = rev(receptor_list.ctrl.ir$`Dendritic cells`$gene),dot.min = 0.05,dot.scale = 5) +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

DotPlot(par.combined.filtered, features = rev(ligand_list.ctrl.ir$`T-cells CD4+CD8+`$gene), dot.min = 0.05,dot.scale = 5) +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(par.combined.filtered, features = rev(receptor_list.ctrl.ir$`T-cells CD4+CD8+`$gene),dot.min = 0.05, dot.scale = 5) +theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


##### create matrix with number of potential pairs
counter = c(0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0)  #create counter vector with 16 slots for 16 cell types
number.of.pairs <- matrix(ncol = 16, nrow = 16, dimnames = list(colnames(Strong.potential.Pairs.ctrl.ir[3:18]), colnames(Strong.potential.Pairs.ctrl.ir[3:18]) )) # create 16x16 matrix to save number of pairs between cell types
for (j in 1:16){
  counter = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  for (i in 1:NROW(Strong.potential.Pairs.ctrl.ir)){
    if(Strong.potential.Pairs.ctrl.ir[i,j+2]>0){
      for(k in 1:16){
        if(Strong.potential.Pairs.ctrl.ir[i,k+21]>0){
          counter[k] <- counter[k] + 1
        }
      }
    }
  }
  number.of.pairs[j,] <- counter[1:16]
}

#### The transposed matrix of 'Number of Pairs' is organized to plot the recipient (receptor) cells in columns and provider (ligands) cells in rows
#### The resulting plot is ligands from X population towards receptors in all other cells
number.of.pairs <- t(number.of.pairs)

#### chord diagrams

write.csv(number.of.pairs, file = "PG Ctrl_IR Number of potential pairs (immune detailed).csv")

# B-cells and M0
border_mat = matrix("black", nrow = 2, ncol = ncol(number.of.pairs))
rownames(border_mat) = c("B-cells", "M0")
colnames(border_mat) = colnames(number.of.pairs)
grid.col = celltype.colors

col_mat = rand_color(length(number.of.pairs), transparency = 0.5)
dim(col_mat) = dim(number.of.pairs)  # to make sure it is a matrix

dev.off()
circos.par(gap.degree=5, gap.after=5)
chordDiagram(number.of.pairs[c(8,14),], grid.col = grid.col, link.lwd = 1, link.lty = 1, link.border = border_mat, 
             symmetric = F, directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.width = 0.1, 
             link.arr.length = 0.1, link.arr.type = "big.arrow",
             annotationTrack = "grid", link.largest.ontop = T, link.arr.lty = 1,grid.border = 1,
             row.col = c("#000000", "olivedrab1"))
circos.clear()

## T CELLS
border_mat = matrix("black", nrow = 5, ncol = ncol(number.of.pairs))
rownames(border_mat) = colnames(number.of.pairs)[9:13]
colnames(border_mat) = colnames(number.of.pairs)

circos.par(gap.degree=5, gap.after=5)
chordDiagram(number.of.pairs[c(9:13),], grid.col = grid.col, link.lwd = 1, link.lty = 1, link.border = border_mat, 
             symmetric = F, directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.width = 0.1, 
             link.arr.length = 0.1, link.arr.type = "big.arrow",
             annotationTrack = "grid", link.largest.ontop = T, link.arr.lty = 1,grid.border = 1,
             row.col = celltype.colors[9:13])
circos.clear()

# DCs and NK cells
border_mat = matrix("black", nrow = 2, ncol = ncol(number.of.pairs))
rownames(border_mat) = colnames(number.of.pairs)[15:16]
colnames(border_mat) = colnames(number.of.pairs)

circos.par(gap.degree=5, gap.after=5)
chordDiagram(number.of.pairs[c(15:16),], grid.col = grid.col, link.lwd = 1, link.lty = 1, link.border = border_mat, 
             symmetric = F, directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.width = 0.1, 
             link.arr.length = 0.1, link.arr.type = "big.arrow",
             annotationTrack = "grid", link.largest.ontop = T, link.arr.lty = 1,grid.border = 1,
             row.col = celltype.colors[15:16])

circos.clear()


### Plot ligands from other cells to immune populations ###
### B-cells and M0
test <- t(number.of.pairs)
border_mat = matrix("black", ncol = ncol(number.of.pairs), nrow = 2)
colnames(border_mat) = colnames(number.of.pairs)
rownames(border_mat) =  c("B-cells",  "M0")
grid.col = celltype.colors

col_mat = rand_color(length(number.of.pairs), transparency = 0.5)
dim(col_mat) = dim(test)  # to make sure it is a matrix

col_mat[c(1:16),c(1:16)]<- '#00000010'
col_mat[8,] <- celltype.colors
col_mat[14,] <- celltype.colors

dev.off()
circos.par(gap.degree=5, gap.after=5)
chordDiagram(test[c(8,14),], grid.col = grid.col, link.lwd = 1, link.lty = 1, link.border = border_mat, 
             symmetric = F, directional = -1, direction.type = c("diffHeight", "arrows"), link.arr.width = 0.1, 
             link.arr.length = 0.1, link.arr.type = "big.arrow",
             annotationTrack = "grid", link.largest.ontop = T, link.arr.lty = 1,grid.border = 1,
             col = col_mat[c(8,14),])
circos.clear()


### T-cells
border_mat = matrix("black", ncol = ncol(number.of.pairs), nrow = 5)
colnames(border_mat) = colnames(number.of.pairs)
rownames(border_mat) =  rownames(test)[9:13]

col_mat[c(1:16),c(1:16)]<- '#00000010'

col_mat[9,] <- celltype.colors
col_mat[10,] <- celltype.colors
col_mat[11,] <- celltype.colors
col_mat[12,] <- celltype.colors
col_mat[13,] <- celltype.colors


dev.off()
circos.par(gap.degree=5, gap.after=5)
chordDiagram(test[c(9:13),], grid.col = grid.col, link.lwd = 1, link.lty = 1, link.border = border_mat, 
             symmetric = F, directional = -1, direction.type = c("diffHeight", "arrows"), link.arr.width = 0.1, 
             link.arr.length = 0.1, link.arr.type = "big.arrow",
             annotationTrack = "grid", link.largest.ontop = T, link.arr.lty = 1,grid.border = 1,
             col = col_mat[c(9:13),])
circos.clear()


### DCs and NKcells
border_mat = matrix("black", ncol = ncol(number.of.pairs), nrow = 2)
colnames(border_mat) = colnames(number.of.pairs)
rownames(border_mat) =  rownames(test)[15:16]

col_mat[c(1:16),c(1:16)]<- '#00000010'

col_mat[15,] <- celltype.colors
col_mat[16,] <- celltype.colors

dev.off()
circos.par(gap.degree=5, gap.after=5)
chordDiagram(test[c(15:16),], grid.col = grid.col, link.lwd = 1, link.lty = 1, link.border = border_mat, 
             symmetric = F, directional = -1, direction.type = c("diffHeight", "arrows"), link.arr.width = 0.1, 
             link.arr.length = 0.1, link.arr.type = "big.arrow",
             annotationTrack = "grid", link.largest.ontop = T, link.arr.lty = 1,grid.border = 1,
             col = col_mat[c(15:16),])
circos.clear()


