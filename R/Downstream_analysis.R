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

### Extracting UnNormalized counts for each cluster and prepare file for ballon/dot Plot##
f1 <-read.table("Gene_list.txt", row.names = 1)
cells.use <- WhichCells(object = test , ident = "0 0")
expr <- GetAssayData(object = test, assay.type = "RNA", slot = "counts")[, cells.use]
expr <- as(Class = 'matrix', object = expr)
dim(expr)
18377   668
temp<- expr
temp<-cbind(temp, Total = rowSums(temp))
mer<- merge(f1,temp, by=0)
Ballon_data<- mer[, c("Row.names","Total")]
Ballon_data$Gonadotrope_Intermediate<-mer[, c("Row.names","Total")]
Ballon_data$Row.names<-NULL
Ballon_data$Total<-NULL


cells.use <- WhichCells(object = test , ident = "0 1")
expr <- GetAssayData(object = test, assay.type = "RNA", slot = "counts")[, cells.use]
expr <- as(Class = 'matrix', object = expr)
dim(expr)
18377   416
temp<- expr
temp<-cbind(temp, Total = rowSums(temp))
mer<- merge(f1,temp, by=0)
Ballon_data$Gonadotrope_lhb<-mer[, c("Row.names","Total")]

cells.use <- WhichCells(object = test , ident = "1")
expr <- GetAssayData(object = test, assay.type = "RNA", slot = "counts")[, cells.use]
expr <- as(Class = 'matrix', object = expr)
dim(expr)
18377   364
temp<- expr
temp<-cbind(temp, Total = rowSums(temp))
mer<- merge(f1,temp, by=0)
Ballon_data$Blood<-mer[, c("Row.names","Total")]


cells.use <- WhichCells(object = test , ident = "2")
expr <- GetAssayData(object = test, assay.type = "RNA", slot = "counts")[, cells.use]
expr <- as(Class = 'matrix', object = expr)
dim(expr)
18377   299
temp<- expr
temp<-cbind(temp, Total = rowSums(temp))
mer<- merge(f1,temp, by=0)
Ballon_data$Stem_subcluster<-mer[, c("Row.names","Total")]



cells.use <- WhichCells(object = test , ident = "3")
expr <- GetAssayData(object = test, assay.type = "RNA", slot = "counts")[, cells.use]
expr <- as(Class = 'matrix', object = expr)
dim(expr)
18377   197
temp<- expr
temp<-cbind(temp, Total = rowSums(temp))
mer<- merge(f1,temp, by=0)
Ballon_data$Stem<-mer[, c("Row.names","Total")]


cells.use <- WhichCells(object = test , ident = "4 0")
expr <- GetAssayData(object = test, assay.type = "RNA", slot = "counts")[, cells.use]
expr <- as(Class = 'matrix', object = expr)
dim(expr)
18377    84
temp<- expr
temp<-cbind(temp, Total = rowSums(temp))
mer<- merge(f1,temp, by=0)
Ballon_data$TSHB<-mer[, c("Row.names","Total")]

cells.use <- WhichCells(object = test , ident = "4 1")
expr <- GetAssayData(object = test, assay.type = "RNA", slot = "counts")[, cells.use]
expr <- as(Class = 'matrix', object = expr)
dim(expr)
18377    53
temp<- expr
temp<-cbind(temp, Total = rowSums(temp))
mer<- merge(f1,temp, by=0)
Ballon_data$GH<-mer[, c("Row.names","Total")]

cells.use <- WhichCells(object = test , ident = "4 2")
expr <- GetAssayData(object = test, assay.type = "RNA", slot = "counts")[, cells.use]
expr <- as(Class = 'matrix', object = expr)
dim(expr)
18377    15
temp<- expr
temp<-cbind(temp, Total = rowSums(temp))
mer<- merge(f1,temp, by=0)
Ballon_data$SMTLA<-mer[, c("Row.names","Total")]

cells.use <- WhichCells(object = test , ident = "5")
expr <- GetAssayData(object = test, assay.type = "RNA", slot = "counts")[, cells.use]
expr <- as(Class = 'matrix', object = expr)
dim(expr)
18377   136
temp<- expr
temp<-cbind(temp, Total = rowSums(temp))
mer<- merge(f1,temp, by=0)
Ballon_data$Fshb<-mer[, c("Row.names","Total")]

cells.use <- WhichCells(object = test , ident = "6")
expr <- GetAssayData(object = test, assay.type = "RNA", slot = "counts")[, cells.use]
expr <- as(Class = 'matrix', object = expr)
dim(expr)
18377   123
temp<- expr
temp<-cbind(temp, Total = rowSums(temp))
mer<- merge(f1,temp, by=0)
Ballon_data$Lacto<-mer[, c("Row.names","Total")]

cells.use <- WhichCells(object = test , ident = "7")
expr <- GetAssayData(object = test, assay.type = "RNA", slot = "counts")[, cells.use]
expr <- as(Class = 'matrix', object = expr)
dim(expr)
18377   115
temp<- expr
temp<-cbind(temp, Total = rowSums(temp))
mer<- merge(f1,temp, by=0)
Ballon_data$Non_endo<-mer[, c("Row.names","Total")]


cells.use <- WhichCells(object = test , ident = "8 0")
expr <- GetAssayData(object = test, assay.type = "RNA", slot = "counts")[, cells.use]
expr <- as(Class = 'matrix', object = expr)
dim(expr)
18377   54
temp<- expr
temp<-cbind(temp, Total = rowSums(temp))
mer<- merge(f1,temp, by=0)
Ballon_data$Melano<-mer[, c("Row.names","Total")]

cells.use <- WhichCells(object = test , ident = "8 1")
expr <- GetAssayData(object = test, assay.type = "RNA", slot = "counts")[, cells.use]
expr <- as(Class = 'matrix', object = expr)
dim(expr)
18377   24
temp<- expr
temp<-cbind(temp, Total = rowSums(temp))
mer<- merge(f1,temp, by=0)
Ballon_data$Cortoco<-mer[, c("Row.names","Total")]

cells.use <- WhichCells(object = test , ident = "9")
expr <- GetAssayData(object = test, assay.type = "RNA", slot = "counts")[, cells.use]
expr <- as(Class = 'matrix', object = expr)
dim(expr)
18377   66
temp<- expr
temp<-cbind(temp, Total = rowSums(temp))
mer<- merge(f1,temp, by=0)
Ballon_data$Lactotrop_sub<-mer[, c("Row.names","Total")]

cells.use <- WhichCells(object = test , ident = "10")
expr <- GetAssayData(object = test, assay.type = "RNA", slot = "counts")[, cells.use]
expr <- as(Class = 'matrix', object = expr)
dim(expr)
18377   15
temp<- expr
temp<-cbind(temp, Total = rowSums(temp))
mer<- merge(f1,temp, by=0)
Ballon_data$Novel<-mer[, c("Row.names","Total")]
write.table(Ballon_data, "Female_Ballon_UnArranged_data3AUG.txt", sep = "\t", quote = F)

head(Ballon_data)
Gonadotrope_lhb.Row.names Gonadotrope_lhb.Total Gonadotrope_Intermediate.Row.names Gonadotrope_Intermediate.Total    Blood.Row.names Blood.Total Stem_subcluster.Row.names Stem_subcluster.Total     Stem.Row.names Stem.Total
1        ENSORLG00000003553                 48244                 ENSORLG00000003553                          89413 ENSORLG00000003553         394        ENSORLG00000003553                 10635 ENSORLG00000003553        725
2        ENSORLG00000013460                   567                 ENSORLG00000013460                            443 ENSORLG00000013460          88        ENSORLG00000013460                   345 ENSORLG00000013460        190
3        ENSORLG00000016928                  3994                 ENSORLG00000016928                           3447 ENSORLG00000016928         618        ENSORLG00000016928                  6048 ENSORLG00000016928       4167
4        ENSORLG00000019556                  1360                 ENSORLG00000019556                           1148 ENSORLG00000019556         155        ENSORLG00000019556                   695 ENSORLG00000019556       1458
5        ENSORLG00000025908                  3805                 ENSORLG00000025908                           3247 ENSORLG00000025908         634        ENSORLG00000025908                 12065 ENSORLG00000025908       1202
6        ENSORLG00000029237                 12581                 ENSORLG00000029237                           9930 ENSORLG00000029237         186        ENSORLG00000029237                  9043 ENSORLG00000029237        605
#Arrange the file and replace total with counts name for ballon/dot plot

#########################################Making Ballon/Dot plot####################################

#read file
balloon<-read.delim("Final_Female_", sep = "\t")
dim(balloon)
105   4
head(balloon)
Row.names Gene_names                     Clusters counts
1 ENSORLG00000003553        lhb Gonadotrope lhb cell cluster  48244
2 ENSORLG00000013460      smtla Gonadotrope lhb cell cluster    567
3 ENSORLG00000016928        prl Gonadotrope lhb cell cluster   3994
4 ENSORLG00000019556        gh1 Gonadotrope lhb cell cluster   1360
5 ENSORLG00000025908       pomc Gonadotrope lhb cell cluster   3805
6 ENSORLG00000029237       fshb Gonadotrope lhb cell cluster  12581
library(ggplot2)
level_cluster <- c("Gonadotrope lhb cell cluster",  "Gonadotrope intermediate cluster", "Gonadotrope fshb cell cluster", "Lactotrope cluster","Lactotrope sub-cluster", "Stem cell cluster","Stem cell sub-cluster", "Somatotrope cluster","Thyrotropes cluster", "Somatolactotrope cluster",    "Corticotrope cluster", "Melanotrope cluster","Non-endocrine cluster", "Blood cell cluster","Novel cluster")
level_gene<-c("lhb", "fshb","tshb", "prl", "pomc","gh1", "smtla")
##Extract color for each cluster
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
color_list <- ggplotColours(n=15)
color_list #check color on this https://www.colorhexa.com/
"#F8766D" "#E58700" "#C99800" "#A3A500" "#6BB100" "#00BA38" "#00BF7D" "#00C0AF" "#00BCD8" "#00B0F6" "#619CFF" "#B983FF" "#E76BF3" "#FD61D1" "#FF67A4

#Levels: Blood cluster (#C99800),  Corticotrop cluster (#E76BF3),  Gonadotrope fshb cell cluster (#00BCD8),  Gonadotrope intermediate cluster (#F8766D),  Gonadotrope lhb cell cluster (#E58700), Lactotrope cluster (#00B0F6), Lactotrope sub-cluster (#FD61D1), Melanotrope cluster (#B983FF), Non-endocrine cluster (#619CFF), Novel cluster (#FF67A4"), Somatolactotrope cluster (#00C0AF) , Somatotrope cluster (#00BF7D), Stem cell cluster (#6BB100), Stem cell sub-cluster (#A3A500), Thyrotropes cluster ((#00BA38))

col_val=c( "#C99800","#E76BF3","#00BCD8","#F8766D" ,"#E58700","#00B0F6", "#FD61D1","#B983FF","#619CFF","#FF67A4","#00C0AF", "#00BF7D", "#6BB100","#A3A500","#00BA38")

pdf("Female_Baloon_3AUG.pdf", width=10, height=10)
p6 <- ggplot(balloon, aes(x = factor(Clusters,level = level_cluster), y =  factor(Gene_names, level = level_gene), size = counts, fill= Clusters)) + scale_fill_manual(values = col_val) + geom_point(shape = 21) + theme_bw() + theme() + scale_size_area(max_size=18) + ggtitle("Marker genes") + labs(x = "Cell clusters", y = "Genes") + theme(axis.text.x = element_text(color = "black", size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain"),axis.text.y = element_text(color = "black", size = 15, angle = 0, hjust = 1, vjust = 0, face = "italic"),axis.title.x = element_text(color = "black", size = 18, angle = 0, hjust = .5, vjust = 0, face = "bold"),axis.title.y = element_text(color = "black", size = 18, angle = 90, hjust = .5, vjust = .5, face = "bold"))+ theme(plot.title = element_text(color = "black", size = 20, hjust = .5,face = "bold"))
p6
dev.off()

************************************************************Figure 1 code Female data#####################
setwd("/net/fs-1/home01/sikh/14_Zoomin_clustering_Seurat")
load("/net/fs-1/home01/sikh/14_Zoomin_clustering_Seurat/sample1_28112019.RData")
dim(qcdata)
4945   28
dim(selected)
2644   31

#Fig1B: Panel1
pdf("Fig1b_Panel1_Female_Aug6.pdf", height = 8, width = 9)
plot(log10(qcdata$req_umis_per_cb), log10(qcdata$mito_fraction), pch = 20, col = "#00000022")
points(log10(selected$req_umis_per_cb), log10(selected$mito_fraction), pch = 20, col = "#00ff0055")
points(log10(selected$req_umis_per_cb)[selected$doublet=="Doublet"], log10(selected$mito_fraction)[selected$doublet=="Doublet"], pch = 6, col = "red")
dev.off()
pdf("Fig1b__Panel2_Female_Aug6.pdf", height = 8, width = 9)
plot(log10(qcdata$req_umis_per_cb),log10(rank(-qcdata$req_umis_per_cb)))
points(log10(qcdata$req_umis_per_cb[qcdata$req_umis_per_cb>1500]),log10(rank(-qcdata$req_umis_per_cb[qcdata$req_umis_per_cb>1500])), col = "red")
dev.off()

selected_blood <- qcdata[qcdata$UMAP_1 >5,]
pdf("Fig1c_Female_Aug6.pdf", height = 8, width = 9)
plot(log10(qcdata$req_umis_per_cb), log10(qcdata$gb), pch = 20, col = "#00000022")
points(log10(selected_blood$req_umis_per_cb), log10(selected_blood$gb), pch = 20, col = "red")
dev.off()
