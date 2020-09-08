# @author <Khadeeja Siddique> siddique.khadeeja@nmbu.no
# This code is used to perform initial clustering, sub-clustering.

library(Seurat)
library(dplyr)

### Initial clustering
Re-Seurat for initial clustering:
S1.finalcm <- CreateSeuratObject(counts = final)
S1.finalcm <- NormalizeData(S1.finalcm, normalization.method = "LogNormalize", scale.factor = 10000)
S1.finalcm<- FindVariableFeatures(S1.finalcm, selection.method = "vst", nfeatures = 2000)
LabelPoints(plot = VariableFeaturePlot(scm), points = head(VariableFeatures(S1.finalcm), 10), repel = T)
S1.finalcm <- ScaleData(S1.finalcm, features = rownames(S1.finalcm))
S1.finalcm <- RunPCA(S1.finalcm, features = VariableFeatures(object = finalcm)) ElbowPlot(finalcm)
S1.finalcm <- FindNeighbors(S1.finalcm, dims = 1:10)
S1.finalcm <- RunUMAP(S1.finalcm, dims = 1:10)
S1.finalcm <- FindClusters(S1.finalcm, resolution = 0.5) #11 clusters
DimPlot(S1.finalcm, label = T)

### Sub-clustering

### Zoomed-in POMC expressing cluster
S1.Pomc <- subset(S1.finalcm, idents = 8)
S1.Pomc<- FindNeighbors(S1.Pomc, dims = 1:10)
S1.Pomc<- FindClusters(S1.Pomc, resolution = 0.2) #2 subclusters
FeaturePlot(S1.Pomc, features = c( "ENSORLG00000016002", "ENSORLG00000005993", "ENSORLG00000006472", "ENSORLG00000012904", "ENSORLG00000002228", "ENSORLG00000025908"))
S1.markers <- FindMarkers(S1.Pomc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #identified DE genes for 2 sub-clusters
S1.top10 <- S1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf("heatmap_female_pomc.pdf", height = 10, width = 10)
DoHeatmap(S1.Pomc, features = S1.top10$gene) 
dev.off()
S1.top5 <- S1.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

pdf("heatmap_female_pomc_TOP5.pdf", height = 10, width = 10)
DoHeatmap(S1.Pomc, features = S1.top5$gene) + NoLegend()
dev.off()

S1.finalcm$sub_cluster <- as.character(Idents(S1.finalcm))
S1.finalcm$sub_cluster[Cells(S1.Pomc)] <- paste("8",Idents(S1.Pomc))

DimPlot(S1.finalcm, group.by = "sub_cluster", label = T)

### Zoomed-in GH expressing cluster
S1.GH <- subset(S1.finalcm, idents = 4)
S1.GH<- FindNeighbors(S1.GH, dims = 1:10)
S1.GH<- FindClusters(S1.GH, resolution = 0.5)   # 3 sub-clusters
FeaturePlot(S1.GH, features = c("ENSORLG00000019556", "ENSORLG00000013460", "ENSORLG00000029251"))
S1.finalcm$sub_cluster[Cells(S1.GH)] <- paste("4",Idents(S1.GH))

S1.markers <- FindMarkers(S1.GH, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #identified DE genes for 3 sub-clusters
S1.top10 <- S1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf("heatmap_female_GH.pdf", height = 10, width = 10)
DoHeatmap(S1.GH, features = S1.top10$gene) 
dev.off()
S1.top5 <- S1.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

pdf("heatmap_female_GH_TOP5.pdf", height = 10, width = 10)
DoHeatmap(S1.GH, features = S1.top5$gene) + NoLegend()
dev.off()

###Zooming Gonadotropes##
S1.Gon <- subset(S1.finalcm, idents = c(0))
S1.Gon <- FindNeighbors(S1.Gon, dims = 1:10)
FeaturePlot(S1.Gon, features = c( "GFP-cassette",  "ENSORLG00000003553", "ENSORLG00000029237", "ENSORLG00000029251" ,"ENSORLG00000001780"))
S1.Gon <- FindClusters(S1.Gon, resolution = 0.3) #2 sub-clusters
DimPlot(S1.Gon)

S1.markers <- FindMarkers(S1.Gon, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #identified DE genes for 2 sub-clusters
S1.top10 <- S1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf("heatmap_female_Gon.pdf", height = 10, width = 10)
DoHeatmap(S1.Gon, features = S1.top10$gene) 
dev.off()
S1.top5 <- S1.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

pdf("heatmap_female_GH_TOP5.pdf", height = 10, width = 10)
DoHeatmap(S1.Gon, features = S1.top5$gene) + NoLegend()
dev.off()

pdf("Female_15clustering.pdf")
test$sub_cluster[Cells(S1.Gon)] <- paste("0",Idents(S1.Gon))
DimPlot(S1.finalcm, group.by = "sub_cluster", label = T) #15 clusters
dev.off()

#Adding to seurat object
S1.finalcm$seurat_clusters<-S1.finalcm$sub_cluster
pident=as.factor(S1.finalcm$seurat_clusters) #15 Levels: 0 0 0 1 1 10 2 3 4 0 4 1 4 2 5 6 7 8 0 8 1 9
S1.finalcm@active.ident=pident #changed it with new 15 clusters
S1.finalcm@active.ident  # Levels: 0 0 0 1 1 10 2 3 4 0 4 1 4 2 5 6 7 8 0 8 1 9
levels(test) <- c("0 0","0 1","1","2","3","4 0","4 1","4 2","5","6","7","8 0","8 1","9","10") #make a order
levels(test)


#########################Manually changes cell identity### Manually change the cell identity because many cells should be in other cluster. In cluster thyrotropes "4 0"  many cells are part of somatolactin so I manually change the cells idenity
Cells_smtla<- c("ACTTTCATCTGGAGCC","TGCGGGTAGATCCGAG","TGCGGGTAGATCCGAG","CATCAGACAATGGTCT","CAGCATACAATCTACG") # 5 cells are part of somatolactotropes
SMTLA.old <- subset(S1.finalcm, idents = "4 2")
dim(SMTLA.old)
S1.finalcm <- SetIdent(object = S1.finalcm, cells = Cells_smtla, value = '4 2')
SMTLA.new <- subset(S1.finalcm, idents = "4 2")
dim(SMTLA.new)

#There are few cells that are part of Thyrotropes## I changed it manually

final.SMTLA.cells<-WhichCells(object = SMTLA.new, expression = ENSORLG00000013460 >= 4)
final.Thyro.cells<-WhichCells(object = SMTLA.new, expression = ENSORLG00000013460 <= 3)
test <- SetIdent(object = test, cells = final.SMTLA.cells, value = '4 2')
test <- SetIdent(object = test, cells = final.Thyro.cells, value = '4 0')
Final.SMTLA.Cluster <- subset(test, idents = "4 2")
dim(Final.SMTLA.Cluster)
18377    15
Final.Thyro.Cluster <- subset(test, idents = "4 0")
dim(Final.Thyro.Cluster)
18377    88

#Changes in somatotroph "4 1" cluster. There are few cells which are part of somatotrop but present in cluster "4 0" thyrotroph
Somato_cluster95cells <- subset(test, idents = "4 1")
dim(Somato_cluster95cells)
18377    42
Tshb_cluster95cells <- subset(test, idents = "4 0")
final.Somato.cells<-WhichCells(object = Tshb_cluster95cells, expression = ENSORLG00000019556 > "50")
test <- SetIdent(object = test, cells = final.Somato.cells, value = '4 1')
Final.Somato.Cluster <- subset(test, idents = "4 1")
dim(Final.Somato.Cluster)
18377    53
Final.Thyro.Cluster <- subset(test, idents = "4 0")
dim(Final.Thyro.Cluster)
18377    84


pdf("Final_cluster15_noLabel.pdf", height = 10, width = 13)
DimPlot(test, reduction = "umap", label = F, pt.size = 0.9)
dev.off()

pdf("Final_cluster15_Label.pdf", height = 10, width = 13)
DimPlot(test, reduction = "umap", label = T, pt.size = 0.9)
dev.off()

#############################Heatmap###

S1.markers <- FindAllMarkers(test, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #identified DE genes for 15 clusters
S1.top10 <- S1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
pdf("heatmap_female_15clusters3AUG.pdf", height = 30, width = 30)
DoHeatmap(test, features = S1.top10$gene) + NoLegend()
dev.off()
write.table(S1.top10, file="cluster15_Top10_marker_genes3AUG.txt", sep = "\t")
S1.top5 <- S1.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

pdf("heatmap_female_15clusters_top5J3AUG.pdf", height = 30, width = 30)
DoHeatmap(test, features = S1.top5$gene) + NoLegend()
dev.off()

saveRDS(test, file = "finalS1.rds")
