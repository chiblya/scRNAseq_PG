library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(Hmisc)
library(circlize)
library(Matrix)
library(viridis)

#######################################
######### MERGE CTRL AND IR ###########
#This script continues from PG IR analysis
#######################################

#Merge and recluster
par.combined <- merge(par.control.filtered, y = par.ir.filtered, add.cell.ids = c("CTRL", "IR"), project = "UTvsIR - PG")
par.combined <- NormalizeData(par.combined)
par.combined <- FindVariableFeatures(par.combined, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
all.genes <- rownames(par.combined)
par.combined <- ScaleData(par.combined, features = all.genes)
par.combined <- RunPCA(par.combined, npcs = 30, verbose = FALSE)
ElbowPlot(par.combined)

par.combined <- FindNeighbors(par.combined, reduction = "pca", dims = 1:12)
par.combined <- FindClusters(par.combined, resolution = 0.8)
par.combined <- RunTSNE(par.combined, reduction = "pca", dims = 1:12)

DimPlot(par.combined, group.by = "seurat_clusters", repel = F,pt.size = 1,label = T,label.size = 6) + theme(axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 16, color = "black"))+NoLegend()
DimPlot(par.combined, group.by = "Treatment", pt.size = 0.25,label = F,label.size = 8, cols = c("#000099", "goldenrod2"))+ theme(axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 16, color = "black"))
DimPlot(par.combined, group.by = "CellType.fixed", pt.size = 1, cols = celltype.colors)+ theme(axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 16, color = "black"))


#### Separate epithelial and immune clusters

DimPlot(par.combined, pt.size = 1,label = F, group.by = "broad.annotation", cols = "Set2") + theme(axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 16, color = "black"))

##### Separate Immune populations to achieve better clustering######
#### Immune clusters
par.combined <- SetIdent(par.combined, value = "broad.annotation")
DimPlot(par.combined)
par.combined$backupID <- par.combined$CellType.fixed

par.immune <- subset(par.combined, idents = "Immune")
par.immune <- NormalizeData(par.immune)
par.immune <- FindVariableFeatures(par.immune, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
all.genes <- rownames(par.immune)
par.immune <- ScaleData(par.immune, features = all.genes)
par.immune <- RunPCA(par.immune, npcs = 30, verbose = FALSE)
ElbowPlot(par.immune)
par.immune <- FindNeighbors(par.immune, reduction = "pca", dims = 1:10)
par.immune <- FindClusters(par.immune, resolution = 0.8)
par.immune <- RunTSNE(par.immune, reduction = "pca", dims = 1:10)

DimPlot(par.immune, group.by = "seurat_clusters",pt.size = 1,label = T,label.size = 6) + theme(axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 16, color = "black"))
DimPlot(par.immune, group.by = "Treatment", pt.size = 0.5,label = F,label.size = 8, cols = c("#000099", "goldenrod2"))+ theme(axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 16, color = "black"))

par.immune <- SetIdent(par.immune, value = "CellType.fixed")
DimPlot(subset(par.immune, idents = c("Acinar 1", "Acinar 2", "Intercalated duct",
                                      "Stromal", "Striated duct"), invert=T), group.by = "CellType.fixed", pt.size = 0.8, cols = "Dark2")+ theme(axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 16, color = "black"))


immune.par.markers <- FindAllMarkers(par.immune, logfc.threshold = 0.3, min.pct = 0.25, only.pos = T)
immune.par.markers <- immune.par.markers[immune.par.markers$p_val_adj<0.05, ]
immune.par.markers.top <- immune.par.markers%>% group_by(cluster) %>% top_n(10, avg_logFC)

##### Separate Epithelial populations to achieve better clustering######
#### Epithelial clusters
par.combined <- SetIdent(par.combined, value = "CellType.fixed")
DimPlot(par.combined)
par.epithelium <- subset(par.combined, idents = c("Acinar 1", "Acinar 2", "Intercalated duct", "Striated duct", "Myoepithelial"))
par.epithelium <- NormalizeData(par.epithelium, normalization.method = "LogNormalize", scale.factor = 10000)
par.epithelium <- FindVariableFeatures(par.epithelium, selection.method = "vst")
all.genes <- rownames(par.epithelium)
par.epithelium <- ScaleData(par.epithelium, features = all.genes)
par.epithelium <- RunPCA(par.epithelium, npcs = 30, verbose = FALSE)
ElbowPlot(par.epithelium)
par.epithelium <- FindNeighbors(par.epithelium, reduction = "pca", dims = 1:5, force.recalc = T)
par.epithelium <- FindClusters(par.epithelium, resolution = 0.5)
par.epithelium <- RunTSNE(par.epithelium, reduction = "pca", dims = 1:5)

DimPlot(par.epithelium, group.by = "seurat_clusters",pt.size = 1,label = T,label.size = 6) + theme(axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 16, color = "black")) + NoLegend()
DimPlot(par.epithelium, group.by = "Treatment", pt.size = 1,label = F,label.size = 8, cols = c("#000099", "goldenrod2"))+ theme(axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 16, color = "black"))
DimPlot(par.epithelium, group.by = "CellType.fixed", pt.size = 1,label = F,label.size = 8, cols = "Dark2")+ theme(axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 16, color = "black"))

###################################################
###################################################
### For figure 5
DimPlot(par.combined, pt.size = 1,cols = celltype.colors)+ theme(axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 16, color = "black"))
DimPlot(par.combined, group.by = "Treatment", pt.size = 0.5,label = F,label.size = 8, cols = c("#000099", "goldenrod2"))+ theme(axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 16, color = "black"))

par.combined$CellType.Tx <- paste0(par.combined$CellType.fixed, "_", par.combined$Treatment)
finalcellnumbers <- as.data.frame(table(par.combined$CellType.Tx))

saveRDS(par.combined, file = "Parotid - Ctrl vs IR - annotated (filtered cells).rds")

################################################################
#### Determine markers that change with IR per cell type #######
################################################################
par.combined<-SetIdent(par.combined, value = "Treatment")
ctrl.vs.ir <- FindMarkers(par.combined, ident.1 = "CTRL", ident.2 = "IR")
ctrl.vs.ir <- ctrl.vs.ir[ctrl.vs.ir$p_val_adj<0.05, ]
ctrl.vs.ir <- ctrl.vs.ir[order(ctrl.vs.ir$avg_logFC), ]
ctrl.vs.ir$gene <- rownames(ctrl.vs.ir)
write.csv(ctrl.vs.ir, file = "Control vs IR markers (with all cells combined).csv")

DoHeatmap(object = subset(par.combined, downsample=100),size = 3,disp.min = -1, features = ctrl.vs.ir$gene,group.by = "Treatment", group.bar = T, label = T,draw.lines = T)+ scale_fill_viridis_c(option = 'B')


par.combined <- SetIdent(par.combined, value = "CellType.Tx")
DimPlot(par.combined)
Idents(par.combined) <- factor(Idents(par.combined), levels = sort(levels(Idents(par.combined)),decreasing = F))
par.combined$CellType.Tx <- Idents(par.combined)
comparisons.list <- list()
j=0
for (i in seq(from = 1,to =  length(levels(Idents(par.combined))), by = 2)) {
  comparisons.list[[j+1]] <- FindMarkers(par.combined, ident.1 = levels(Idents(par.combined))[i], ident.2 = levels(Idents(par.combined))[i+1], only.pos = F, logfc.threshold = 0.5, min.pct = 0.25)  
  comparisons.list[[j+1]]$gene <- rownames(comparisons.list[[j+1]])
  names(comparisons.list[[j+1]])[3] <- levels(Idents(par.combined))[i]
  names(comparisons.list[[j+1]])[4] <- levels(Idents(par.combined))[i+1]
  comparisons.list[[j+1]]<-comparisons.list[[j+1]][comparisons.list[[j+1]]$p_val_adj <0.05, ]
  comparisons.list[[j+1]]<-comparisons.list[[j+1]][order(comparisons.list[[j+1]]$avg_logFC), ]
  j=j+1
}

#### Extract number of IR-disregulated genes per cell and graph (For figure 6C-D)
number.of.degs.per.cell.post.ir <- c()
test <- plyr::ldply(comparisons.list, rbind)
test <- test[,-c(1:3, 6:7)]
test <- test[,-c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31)]
test2 <- as.data.frame(sapply(test, function(x) sum(!(is.na(x)))))
names(test2)[1] <- "degs"
rownames(test2)<-sub("_IR", "", rownames(test2))
names(comparisons.list) <- rownames(test2) # add list names for later

par(mar=c(9,4,2,2))
barplot(test2$degs,names.arg = rownames(test2),
        ylab="Number of IR-DEGs",
        border="black", col = "blue", las=2, axis.lty = 1)

## create new sorting order (alphabetically to ensure coherence between all exported files)
par.combined <- SetIdent(par.combined, value = "CellType.fixed")
Idents(par.combined) <- factor(Idents(par.combined), levels = sort(levels(Idents(par.combined)),decreasing = F))
DimPlot(par.combined)

## CD4+ and CD4+CD8+ are flipped, so we'll reorder comparisons.list 
## to match the desired sorting order
test.list <- vector(mode = "list", length = 16)
names(test.list) <- levels(Idents(par.combined))
comparisons.list[names(test.list)]
comparisons.list <- comparisons.list[names(test.list)]
## remove temporary file
remove(test.list)

#### Export tables 
## Maintain same order in IDENTS
par.combined <- SetIdent(par.combined, value = "CellType.fixed")
Idents(par.combined) <- factor(Idents(par.combined), levels = sort(levels(Idents(par.combined)),decreasing = F))
DimPlot(par.combined)

for (i in 1:length(comparisons.list)) {
  identity <- levels(Idents(par.combined))[[i]]
  write.csv(x = comparisons.list[[i]],file = paste0("~/Collab - Kirsten/Parotid Manuscript (KHL collab)/PG control process/Gene Lists/PAR_",identity, " ctrl vs ir.csv"),row.names = T)
}

## The following code exports violin plots for the top up and downregulated genes
## per cell type (if less than 10 are present, only those are shown)
plot.grid <- list()
par.combined <- SetIdent(par.combined, value = "CellType.fixed")
Idents(par.combined) <- factor(Idents(par.combined), levels = sort(levels(Idents(par.combined)),decreasing = F))
DimPlot(par.combined)

for (i in 1:16) {
  identity <- levels(Idents(par.combined))[[i]]
  #upregulated genes
  comparisons.list[[i]]<-comparisons.list[[i]][order(comparisons.list[[i]]$avg_logFC), ]
  plot.list <- list()
  if(sum(comparisons.list[[i]]$avg_logFC<0)>0){
    for (j in 1:sum(comparisons.list[[i]]$avg_logFC<0)) {
      plot.list[[j]] <- VlnPlot(par.combined, idents = identity, features = comparisons.list[[i]]$gene[j],pt.size = 0,cols = c("#000099", "#FF6666"),split.by = "Treatment") +NoLegend() +theme(axis.title = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(size = 18,face = "italic")) +NoAxes()
    }
    if(length(plot.list)>10){
      ggsave(filename = paste0(identity," upregulated genes.png"),dpi = 300, units = "in", width = 20,height = 1.5,plot = plot_grid(plotlist = plot.list[1:10],ncol = 20)) 
    }
    if(length(plot.list)<=10){
      ggsave2(filename = paste0(identity," upregulated genes.png"),dpi = 300, units = "in", width = 20,height = 1.5,plot = plot_grid(plotlist = plot.list,ncol = 20)) 
    }
  }  
  #downregulated genes
  comparisons.list[[i]]<-comparisons.list[[i]][order(comparisons.list[[i]]$avg_logFC,decreasing = T), ]
  plot.list <- list()
  if(sum(comparisons.list[[i]]$avg_logFC>0)>0){
    for (j in 1:sum(comparisons.list[[i]]$avg_logFC>0)) {
      plot.list[[j]] <- VlnPlot(par.combined, idents = identity, features = comparisons.list[[i]]$gene[j],pt.size = 0,cols = c("#000099", "#FF6666"),split.by = "Treatment") +NoLegend() +theme(axis.title = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(size = 18,face = "italic")) +NoAxes()
    }
    if(length(plot.list)>10){
      ggsave2(filename = paste0(identity," downregulated genes.png"),dpi = 300, units = "in", width = 20,height = 1.5,plot = plot_grid(plotlist = plot.list[1:10],ncol = 20)) 
    }
    if(length(plot.list)<=10){
      ggsave2(filename = paste0(identity," downregulated genes.png"),dpi = 300, units = "in", width = 20,height = 1.5,plot = plot_grid(plotlist = plot.list,ncol = 20)) 
    }
  }}

# test code to determine how many unique IR-disregulated genes are identified
##genes <- c()
#for (i in 1:length(comparisons.list)) {
#  genes <- unique(c(genes, comparisons.list[[i]]$gene[1:10]))
#}


 
## determine list of genes that are statistically different in specific cell types
## ctrl vs IR
par.combined <- SetIdent(par.combined, value = "CellType.fixed")
Idents(par.combined) <- factor(Idents(par.combined), levels = sort(levels(Idents(par.combined)),decreasing = F))
DimPlot(par.combined)
list.disgenes.unique <- list()

for (i in 1:length(comparisons.list)) {
  list.disgenes.unique[[i]] <- comparisons.list[[i]]$gene
  for (j in 1:length(comparisons.list)) {
    if (j!=i) {
      list.disgenes.unique[[i]] <- list.disgenes.unique[[i]][!(list.disgenes.unique[[i]] %in% comparisons.list[[j]]$gene)]
    }
  }
}
names(list.disgenes.unique) <- levels(Idents(par.combined)) #matches comparisons list
list.disgenes.unique.table <- plyr::ldply(list.disgenes.unique, rbind)


## plot "unique" degs in epithelial clusters
## these genes may also be expressed by multiple cell types but are uniquely disrupted by 
## IR in the shown cluster

DotPlot(par.combined, features = list.disgenes.unique[[1]],scale.max = 100, col.min = 0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(par.combined, features = list.disgenes.unique[[6]],dot.min = 0.1, col.min = 0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(par.combined, features = list.disgenes.unique[[10]],dot.min = 0.1, col.min = 0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

## plot "unique" degs in T-cell clusters
## these genes may also be expressed by multiple cell types but are uniquely disrupted by 
## IR in the shown cluster

DotPlot(par.combined, features = list.disgenes.unique[[12]],scale.max = 100, col.min = 0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(par.combined, features = list.disgenes.unique[[13]],scale.max = 100, col.min = 0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(par.combined, features = list.disgenes.unique[[14]],scale.max = 100, col.min = 0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(par.combined, features = list.disgenes.unique[[15]],scale.max = 100, col.min = 0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(par.combined, features = list.disgenes.unique[[16]],scale.max = 100, col.min = 0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

## plot "unique" degs in Other Immune clusters
## these genes may also be expressed by multiple cell types but are uniquely disrupted by 
## IR in the shown cluster

DotPlot(par.combined, features = list.disgenes.unique[[3]],scale.max = 100, col.min = 0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(par.combined, features = list.disgenes.unique[[4]],scale.max = 100, col.min = 0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
DotPlot(par.combined, features = list.disgenes.unique[[7]],scale.max = 100, col.min = 0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
