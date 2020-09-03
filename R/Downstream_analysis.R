# @author <Khadeeja Siddique> siddique.khadeeja@nmbu.no
# This code is used to perform initial clustering, sub-clustering, to generate dotplots and heatmaps.

setwd("/net/fs-1/home01/sikh/14_Zoomin_clustering_Seurat")
load("/net/fs-1/home01/sikh/14_Zoomin_clustering_Seurat/sample1_28112019.RData")
library(Seurat)
library(dplyr)
S1.finalcm<-finalcm
S1.finalcm <- FindClusters(S1.finalcm, resolution = 0.5) #11 clusters
DimPlot(S1.finalcm, label = T)

### Zooming POMCS###
S1.Pomc <- subset(S1.finalcm, idents = 8)
S1.Pomc<- FindNeighbors(S1.Pomc, dims = 1:10)
S1.Pomc<- FindClusters(S1.Pomc, resolution = 0.2) #2 subclusters
FeaturePlot(S1.Pomc, features = c( "ENSORLG00000016002", "ENSORLG00000005993", "ENSORLG00000006472", "ENSORLG00000012904", "ENSORLG00000002228", "ENSORLG00000025908"))
test= S1.finalcm
test$sub_cluster <- as.character(Idents(test))
test$sub_cluster[Cells(S1.Pomc)] <- paste("8",Idents(S1.Pomc))
DimPlot(test, group.by = "sub_cluster", label = T)

###Zooming GH###
S1.GH <- subset(S1.finalcm, idents = 4)
S1.GH<- FindNeighbors(S1.GH, dims = 1:10)
S1.GH<- FindClusters(S1.GH, resolution = 0.5)   #3sub-clusters
FeaturePlot(S1.GH, features = c("ENSORLG00000019556", "ENSORLG00000013460", "ENSORLG00000029251"))
test$sub_cluster[Cells(S1.GH)] <- paste("4",Idents(S1.GH))
DimPlot(test, group.by = "sub_cluster", label = T)

###Zooming Gonadotropes##
S1.Gon <- subset(S1.finalcm, idents = c(0))
S1.Gon <- FindNeighbors(S1.Gon, dims = 1:10)
FeaturePlot(S1.Gon, features = c( "GFP-cassette",  "ENSORLG00000003553", "ENSORLG00000029237", "ENSORLG00000029251" ,"ENSORLG00000001780"))
S1.Gon <- FindClusters(S1.Gon, resolution = 0.3) #2 sub-clusters
DimPlot(S1.Gon)
pdf("FinalFemale_15clustering.pdf")
test$sub_cluster[Cells(S1.Gon)] <- paste("0",Idents(S1.Gon))
DimPlot(test, group.by = "sub_cluster", label = T) #15 clusters
dev.off()

test$seurat_clusters<-test$sub_cluster
pident=as.factor(test$seurat_clusters) #15 Levels: 0 0 0 1 1 10 2 3 4 0 4 1 4 2 5 6 7 8 0 8 1 9
test@active.ident #11 clusters old: Levels: 0 1 2 3 4 5 6 7 8 9 10
test@active.ident=pident #changed it with new 15 clusters
test@active.ident  # Levels: 0 0 0 1 1 10 2 3 4 0 4 1 4 2 5 6 7 8 0 8 1 9
levels(test)
[1] "0 0" "0 1" "1"   "10"  "2"   "3"   "4 0" "4 1" "4 2" "5"   "6"   "7"   "8 0" "8 1" "9"
levels(test) <- c("0 0","0 1","1","2","3","4 0","4 1","4 2","5","6","7","8 0","8 1","9","10") #make a order
levels(test)
[1] "0 0" "0 1" "1"   "2"   "3"   "4 0" "4 1" "4 2" "5"   "6"   "7"   "8 0" "8 1" "9"   "10"


############################################################Manually changes cell identity###
#Manually change the cell identity because many cells should be in other cluster. In cluster thyrotropes "4 0"  many cells are part of somatolactin so I manually change the cells idenity
Cells_smtla<- c("ACTTTCATCTGGAGCC","TGCGGGTAGATCCGAG","TGCGGGTAGATCCGAG","CATCAGACAATGGTCT","CAGCATACAATCTACG") # 5 cells are part of somatolactotropes
SMTLA.old <- subset(test, idents = "4 2")
dim(SMTLA.old)
18377    19
test <- SetIdent(object = test, cells = Cells_smtla, value = '4 2')
SMTLA.new <- subset(test, idents = "4 2")
dim(SMTLA.new)
18377    23

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

###Fshb###
Fshb.Cluster5 <- subset(test, idents = "5")
final.fshb.cells<-WhichCells(object = Fshb.Cluster5)
expr <- GetAssayData(object = Fshb.Cluster5, assay.type = "RNA", slot = "counts")[, final.fshb.cells]
expr <- as(Class = 'matrix', object = expr)
write.csv(expr, "Fshb_cluster5.csv")

final.fshb.cells<-WhichCells(object = Fshb.Cluster5, expression = ENSORLG00000029237 > 3)
test <- SetIdent(object = test, cells = final.fshb.cells, value = 'New_Fshb')
DimPlot(test, label = T)
test <- RenameIdents(object = test, 'New_Fshb' = '5', '5' = '0 0')
DimPlot(test, label = T)

## Now for Lhb and intermediate
lhb.Cluster01 <- subset(test, idents = "0 1")
lhb.Initial.cells<-WhichCells(object = lhb.Cluster01)
expr <- GetAssayData(object = lhb.Cluster01, assay.type = "RNA", slot = "counts")[, lhb.Initial.cells]
expr <- as(Class = 'matrix', object = expr)
write.csv(expr, "lhb_cluster01.csv")

final.lhb.cells<-WhichCells(object = lhb.Cluster01, expression = ENSORLG00000003553 > 3)
test <- SetIdent(object = test, cells = final.lhb.cells, value = 'New_lhb')
DimPlot(test, label = T)
test <- RenameIdents(object = test, 'New_lhb' = '0 1', '0 1' = '0 0')
DimPlot(test, label = T)

Intermed.Cluster00 <- subset(test, idents = "0 0")
lntermed.Initial.cells<-WhichCells(object = Intermed.Cluster01)
expr <- GetAssayData(object = Intermed.Cluster01, assay.type = "RNA", slot = "counts")[, lntermed.Initial.cells]
expr <- as(Class = 'matrix', object = expr)
write.csv(expr, "lntermediate_cluster00.csv")

Fshb_new_cells<-WhichCells(object = Intermed.Cluster00, expression = ENSORLG00000029237 > "50")
test <- SetIdent(object = test, cells = Fshb_new_cells, value = 'New2_Fshb')

Lhb_new_cells<-WhichCells(object = Intermed.Cluster00, expression =  ENSORLG00000003553 > "50")
test <- SetIdent(object = test, cells = Lhb_new_cells, value = 'New_Lhb')

Tshb_new_cells<-WhichCells(object = Intermed.Cluster00, expression = ENSORLG00000029251 > "40")
test <- SetIdent(object = test, cells = Tshb_new_cells, value = 'New_Tshb')

test <- RenameIdents(object = test, 'New2_Fshb' = '5', 'New_Lhb' = '0 1', 'New_Tshb'='4 0')
DimPlot(test, label = T)

levels(test) <- c("0 0","0 1","1","2","3","4 0","4 1","4 2","5","6","7","8 0","8 1","9","10") #make a order
levels(test)

pdf("Final_cluster15_noLabelAug3.pdf", height = 10, width = 13)
DimPlot(test, reduction = "umap", label = F, pt.size = 0.9)
dev.off()

pdf("Final_cluster15_LabelAug3.pdf", height = 10, width = 13)
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
