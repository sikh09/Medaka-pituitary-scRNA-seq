# @author <Khadeeja Siddique> siddique.khadeeja@nmbu.no
# This code is used to remove the duplets from scRNA seq data and updated quality control step.

Library(DoubletFinder)
scm <- CreateSeuratObject(counts = cm)
scm <- NormalizeData(scm, normalization.method = "LogNormalize", scale.factor = 10000)
scm <- FindVariableFeatures(scm, selection.method = "vst", nfeatures = 2000)
LabelPoints(plot = VariableFeaturePlot(scm), points = head(VariableFeatures(scm), 10), repel = T)
scm <- ScaleData(scm, features = rownames(scm))
scm <- RunPCA(scm, features = VariableFeatures(object = scm)) DimHeatmap(scm, dims = 1:2, cells = 1000, balanced = T) ElbowPlot(scm)
scm <- RunUMAP(scm, dims = 1:15) #10-15 informative PCA dimensions DimPlot(scm, reduction = "umap")
df_sweep <- paramSweep_v3(scm, PCs = 1:15, sct = F) df_sweep_summary <- summarizeSweep(df_sweep, GT = F) find.pK(df_sweep_summary)
scm <- doubletFinder_v3(scm, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = 0.02 * 2644, reuse.pANN = F, sct = F)

# 2% doublet rate for ~2500 cells as per the 10x paper. Reference
doublets <- data.frame(rownames(scm@meta.data), scm@meta.data$DF.classifications_0.25_0.01_52.88, scm@meta.data$pANN_0.25_0.01_52.88)
colnames(doublets) <- c("barcode", "doublet", "pANN")
umap <- data.frame(scm@reductions$umap@cell.embeddings, row.names = rownames(scm@reductions$umap@cell.embeddings))
selected <- merge(selected, doublets, by.x = "barcode", by.y = "barcode") selected <- merge(selected, umap, by.x = "barcode", by.y = 0)
# Update the (selected) QC table with new UMAP & doublets:
write.table(selected, file ="sample1_postqc_25112019.txt", na = "0", row.names = T, sep = "\t")

plot(selected$UMAP_1.y, selected$UMAP_2.y, pch = ".") points(selected$UMAP_1.y[selected$doublet=="Doublet"], selected$UMAP_2.y[selected$doublet=="Doublet"], col="green", pch = 6) # triangles
points(selected$UMAP_1.y[selected$req_umis_per_cb>25000], selected$UMAP_2.y[selected$req_umis_per_cb>25000], col="red", pch = 4) # *
final <- cm[selected$doublet == "Singlet"]
