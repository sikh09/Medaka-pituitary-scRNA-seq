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

S1.finalcm$sub_cluster <- as.character(Idents(S1.finalcm)) # Adding information
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

#Adding the sub-clusters into seurat object
S1.finalcm$seurat_clusters<-S1.finalcm$sub_cluster
pident=as.factor(S1.finalcm$seurat_clusters) #15 Levels: 0 0 0 1 1 10 2 3 4 0 4 1 4 2 5 6 7 8 0 8 1 9
S1.finalcm@active.ident=pident #changed it with new 15 clusters
S1.finalcm@active.ident  # Levels: 0 0 0 1 1 10 2 3 4 0 4 1 4 2 5 6 7 8 0 8 1 9
levels(S1.finalcm) <- c("0 0","0 1","1","2","3","4 0","4 1","4 2","5","6","7","8 0","8 1","9","10") #make a order
S1.data<- S1.finalcm
# Assign the names to cell clusters
new.cluster.ids <- c("11. Red blood cells" ,  "9. Lhb-gonadotropes" , "8. Fshb-gonadotropes", "10. Gonadotrope-like"  ,    "7. Thyrotropes"  ,     "6. Somatotrope" ,
                     "5. Somatolactotrope" , "1. Melanotrope"   ,    "2. Corticotrope"   ,   "12. Macrophages"    ,   "3. Lactotrope"  ,  "4. Lactotrope"   ,
                     "13. Uncharacterized", "14. Uncharacterized","15. cga-expressing cells")
names(new.cluster.ids) <- levels(S1.data)
S1.data <- RenameIdents(S1.data, new.cluster.ids)
color_S1 = c("#B983FF" ,  "orchid4", "lightcyan3" , "#FD61D1"  , "green2" ,"yellow","orange",  "#00BCD8", "peru",  "#F8766D"  ,  "firebrick3" , "rosybrown1", "#A3A500" ,"#6BB100")
p2<- DimPlot(S1.data, reduction = "umap", cols = color_S1,label = F, pt.size = 0.7)
p2<- p2 + theme(legend.position="bottom",  panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(hjust = 0.5)) + labs(title = "Female Clustering")
p2
pdf("Final_female_clustering.pdf", height = 10, width = 12)
p2
dev.off()

######### Heatmap ##########

library(dplyr)
cluster.markers <- FindAllMarkers(S1.data, min.pct = 0.25)
S1.top5 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
write.table(S1.top5, file="01_List_top5HeatmapGen.txt", sep = "\t")
dim(S1.top5) # 75 5


pdf("01_Female_heatmap.pdf", width = 6, height = 8)
DoHeatmap(S1.data, features = S1.top5$gene, cells = WhichCells(S1.data, idents =  c("1. Melanotropes" ,"2. Corticotropes", "3. Lactotropes",  "4. Lactotropes", "5. Somatolactotropes" , "6. Somatotropes" , "7. Thyrotropes","8. Fshb gonadotropes", "9. Lhb gonadotropes" , "10. Gonadotropes"  ,"11. Red blood cells" , "12. Macrophages","13. Uncharacterized", "14. Uncharacterized","15. Uncharacterized")), group.colors= c("#B983FF"   ,   "orchid4"   ,   "lightcyan3"    , "#FD61D1"    , "green2"  ,"yellow","orange"   ,    "#00BCD8", "peru",  "#F8766D"  ,   "firebrick3" , "rosybrown1", "#A3A500" ,"#6BB100", "blue") , label = F) + NoLegend()
dev.off()
# for male data
S2.cluster.markers <- FindAllMarkers(S2.data, min.pct = 0.25)
S2.top5 <- S2.cluster.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
write.table(S2.top5, file="02_List_Male_top5HeatmapGen.txt", sep = "\t")

pdf("02_Male_heatmapt.pdf", width = 6, height = 8)
DoHeatmap(S2.data, features = S2.top5$gene, cells = WhichCells(S2.data, idents =  c("1. Melanotropes" ,"2. Corticotropes", "3. Lactotropes",  "4. Lactotropes", "5. Somatolactotropes" , "6. Somatotropes" , "7. Thyrotropes","8. Fshb gonadotropes", "9. Lhb gonadotropes" , "10. Gonadotropes"  ,"11. Red blood cells" , "12. Macrophages","13. Uncharacterized", "14. Uncharacterized","15. Uncharacterized","16. Uncharacterized" )), group.colors= c("#B983FF"   ,   "orchid4"   ,   "lightcyan3"    , "#FD61D1"    , "green2"  ,"yellow","orange"   ,    "#00BCD8", "peru",  "#F8766D"  ,   "firebrick3" , "rosybrown1", "#A3A500" ,"#6BB100", "blue", "turquoise4") , label = F)+ NoLegend()
dev.off()


############ Ballon Plot ##########
library(ggpubr)
library(ggplot2)

f1 <-read.table("Gene_list.txt", row.names = 1) #contained the list of hormone-producing gene marker 
S1_Avg_Exp<- AverageExpression(S1.data) # taking the normalized avg expression of gene per cluster
mer<- merge(S1_Avg_Exp$RNA, f1, by=0)
write.table(mer, "Female_ballon_data.txt", sep = "\t") #arrabge it into 4 colums
S1.ballon.data<- read.delim("Female_ballon_data.txt", sep = "\t")
level_gene<-c( "cga","lhb", "fshb","tshb","gh1","smtla","prl", "pomc")
col_val<- c("#B983FF"   ,   "#F8766D"  ,"blue",            "orchid4"   ,   "lightcyan3"    , "#FD61D1"    , "green2"  ,"yellow","orange"   ,    "#00BCD8", "peru")
level_cluster_M= c("1. Melanotropes" ,"2. Corticotropes", "3. Lactotropes",  "4. Lactotropes", "5. Somatolactotropes" , "6. Somatotropes" , "7. Thyrotropes","8. Fshb-gonadotropes", "9. Lhb-gonadotropes" , "10. Gonadotrope-like"  , "11. Red blood cells" , "12. Macrophages","13. Uncharacterized", "14. Uncharacterized","15. cga-expressing cells", "16. Uncharacterized")
level_cluster_F= c("1. Melanotropes" ,"2. Corticotropes", "3. Lactotropes",  "4. Lactotropes", "5. Somatolactotropes" , "6. Somatotropes" , "7. Thyrotropes","8. Fshb-gonadotropes", "9. Lhb-gonadotropes" , "10. Gonadotrope-like"  , "11. Red blood cells" , "12. Macrophages","13. Uncharacterized", "14. Uncharacterized","15. cga-expressing cells")
S2.ballon.data<- read.delim("Final_male_ballon_data", sep = "\t")
p3 <- ggplot(S1.ballon.data, aes(x = factor(Clusters,level = level_cluster_F), y = factor(Gene_names, level = level_gene))) + geom_point(data=S1.ballon.data[which(S1.ballon.data$counts <= 184),], aes(size=counts), fill='gray', color='gray', shape=21) + geom_point(data=S1.ballon.data[which(S1.ballon.data$counts > 185),],aes(fill=Clusters, size=counts), color = "white",shape=21)+ scale_fill_manual(values = col_val) + guides(fill=guide_legend(ncol=2)) + theme_bw() + theme() + scale_size_area(max_size=18) + labs(x = "Cell types", y = "Gene expression") + theme(axis.text.x = element_text(color = "black", size = 15, angle = 90, hjust=0.95,vjust=0.2, face = "plain"),axis.text.y = element_text(color = "black", size = 15, angle = 0, hjust = 1, vjust = 0, face = "italic"),axis.title.x = element_text(color = "black", size = 18, angle = 0, hjust = .5, vjust = 0, face = "bold"),axis.title.y = element_text(color = "black", size = 18, angle = 90, hjust = .5, vjust = .5, face = "bold"))+ theme(plot.title = element_text(color = "black", size = 20, hjust = .5,face = "bold"))
p3= p3 + theme(legend.position = "none")+ labs(title = "Female") + theme(
panel.grid.major.x = element_blank()
) # to remove y axis grid
pdf("Female_Ballon_plot.pdf", width = 12, height = 10)
p3
dev.off()
p4 <- ggplot(S2.ballon.data, aes(x = factor(Clusters,level = level_cluster_M), y = factor(Gene_names, level = level_gene))) + geom_point(data=S2.ballon.data[which(S2.ballon.data$counts <= 184),], aes(size=counts), fill='gray', color='gray', shape=21) + geom_point(data=S2.ballon.data[which(S2.ballon.data$counts > 185),],aes(fill=Clusters, size=counts), color = "white", shape=21)+ scale_fill_manual(values = col_val) + guides(fill=guide_legend(ncol=2)) + theme_bw() + theme() + scale_size_area(max_size=18) + labs(x = "Cell types", y = "Gene expression") + theme(axis.text.x = element_text(color = "black", size = 15, angle = 90, hjust=0.95,vjust=0.2,  face = "plain"),axis.text.y = element_text(color = "black", size = 15, angle = 0, hjust = 1, vjust = 0, face = "italic"),axis.title.x = element_text(color = "black", size = 18, angle = 0, hjust = .5, vjust = 0, face = "bold"),axis.title.y = element_text(color = "black", size = 18, angle = 90, hjust = .5, vjust = .5, face = "bold"))+ theme(plot.title = element_text(color = "black", size = 20, hjust = .5,face = "bold"))
p4=p4+ theme(legend.position = "none")+ labs(title = "Male") + theme(
panel.grid.major.x = element_blank()
)
pdf("Male_Ballon_plot.pdf", width = 12, height = 10)
p4
dev.off()

pdf("Final_ballon_plo.pdf", width = 17.5, height = 8)
ggarrange(p3, p4,
#labels = c("A", "B"),
          ncol = 2, nrow = 1)
dev.off()                   
