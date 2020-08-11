# @author <Khadeeja Siddique> siddique.khadeeja@nmbu.no
# This code is used for preliminary visulization and quality control of scRNA seq data

library(dropestr)
library(Seurat)
setwd("F://scmedaka/sample1_11072019/dropestnew")
# Visulisation
data <- readRDS("sample1_25112019.rds")
cm <- as.matrix(data$cm)
scm <- CreateSeuratObject(counts = cm)
scm <- NormalizeData(scm, normalization.method = "LogNormalize", scale.factor = 10000)
scm <- FindVariableFeatures(scm, selection.method = "vst", nfeatures = 2000)
LabelPoints(plot = VariableFeaturePlot(scm), points = head(VariableFeatures(scm), 10), repel = T)
scm <- ScaleData(scm, features = rownames(scm))
scm <- RunPCA(scm, features = VariableFeatures(object = scm))
ElbowPlot(scm)
scm <- RunUMAP(scm, dims = 1:15) #10-15 informative PCA dimensions DimPlot(scm, reduction = "umap")
umap <- data.frame(scm@reductions$umap@cell.embeddings, row.names = rownames(scm@reductions$umap@cell.embeddings))

# Quality Control
readspercell <- cbind(data$mean_reads_per_umi[order(names(data$mean_reads_per_umi))], data$aligned_reads_per_cell[order(names(data$aligned_reads_per_cell))], data$aligned_umis_per_cell[order(names(data$aligned_umis_per_cell))], data$requested_umis_per_cb[order(names(data$requested_umis_per_cb))], data$requested_reads_per_cb[order(names(data$requested_reads_per_cb))])
mito <- data$reads_per_chr_per_cells$Intergenic$MT[order(rownames(data$reads_per_c hr_per_cells$Intergenic))] + data$reads_per_chr_per_cells$Exon$MT[order(rownames(data$reads_per_chr_per _cells$Exon))]
exons <- rowSums(data$reads_per_chr_per_cells$Exon)[order(rownames(data$reads_per_c hr_per_cells$Exon))] - data$reads_per_chr_per_cells$Exon$MT[order(rownames(data$reads_per_chr_per _cells$Exon))]
intergenic <- rowSums(data$reads_per_chr_per_cells$Intergenic)[order(rownames(data$reads _per_chr_per_cells$Intergenic))] - data$reads_per_chr_per_cells$Intergenic$MT[order(rownames(data$reads_per_c hr_per_cells$Intergenic))]
introns <- rowSums(data$reads_per_chr_per_cells$Intron)[order(rownames(data$reads_per _chr_per_cells$Intron))] # no introns in the mitochondria
qcdata <- cbind(mito, exons, intergenic)
qcdata <- merge(qcdata, introns, by.x = 0, by.y = 0, all = T)
qcdata <- merge(qcdata, data.frame(readspercell), by.x = "Row.names", by.y = 0)
qcdata[is.na(qcdata)] <- 0
colnames(qcdata) <- c("barcode", "mito", "exons", "intergenic", "introns", "reads_per_umi", "reads_per_cell", "umis_per_cell", "req_umis_per_cb", "req_reads_per_cb")
qcdata <- cbind(qcdata, qcdata$mito/qcdata$reads_per_cell, qcdata$exons/qcdata$reads_per_cell, qcdata$intergenic/qcdata$reads_per_cell, qcdata$introns/qcdata$reads_per_cell)
colnames(qcdata)[11:14] <- c("mito_fraction", "exon_fraction", "intergenic_fraction", "intron_fraction")
rownames(qcdata) <- qcdata$barcode
genes_per_cell <- data.frame(colnames(data$cm_raw), data$cm_raw@p[2:74339] - data$cm_raw@p[1:74338])
colnames(genes_per_cell) <- c("barcode", "genes_per_cell")
qcdata <- merge(qcdata, genes_per_cell, by.x = "barcode", by.y = "barcode")
qcdata <- cbind(qcdata, qcdata$umis_per_cell/qcdata$genes_per_cell) colnames(qcdata)[16] <- c("umis_per_gene")
qcdata <- merge(qcdata, umap, by.x = "barcode", by.y = 0) tmpexpr <- as.data.frame(t(as.matrix(scm@assays$RNA@data)))
qcexpr <- data.frame(tmpexpr$'GFP-cassette', tmpexpr$ENSORLG00000019556, tmpexpr$ENSORLG00000025908, tmpexpr$ENSORLG00000003553, tmpexpr$ENSORLG00000029237, tmpexpr$ENSORLG00000016928, tmpexpr$ENSORLG00000013460, tmpexpr$ENSORLG00000003046, tmpexpr$ENSORLG00000029251, tmpexpr$ENSORLG00000022598)
colnames(qcexpr) <- c("gfp", "gh", "pomc", "lhb", "fshb", "prl", "smtl", "gb", "tshb", "cga")
rownames(qcexpr) <- colnames(scm@assays$RNA@data)
qcdata <- merge(qcdata, qcexpr, by.x = 0, by.y = 0)
write.table(qcdata, file ="sample1_qc_25112019.txt", na = "0", row.names = T, sep = "\t")

# Final cutoff decisions
plot(density(qcdata$mito_fraction))
hist(qcdata$mito_fraction, breaks = 60)

# Mitochondrial fraction cutoff: <5%
plot(log10(rank(-qcdata$req_umis_per_cb)), log10(qcdata$req_umis_per_cb)) points(log10(rank(-qcdata$req_umis_per_cb[qcdata$req_umis_per_cb>1500])), log10(qcdata$req_umis_per_cb[qcdata$req_umis_per_cb>1500]), col = "red")

# UMIs per cellular barcode cutoff: 1500
plot(qcdata$UMAP_1, qcdata$UMAP_2, pch = 20, col = "#00000011", main = "excluded based on mitochondrial fraction (>5%)") points(qcdata$UMAP_1[qcdata$mito_fraction>0.05], qcdata$UMAP_2[qcdata$mito_fraction>0.05], pch = 20, col = "#ff000022")
plot(qcdata$UMAP_1, qcdata$UMAP_2, pch = 20, col = "#00000011", main = "excluded based on UMI count (<1500)") points(qcdata$UMAP_1[qcdata$req_umis_per_cb<1500], qcdata$UMAP_2[qcdata$req_umis_per_cb<1500], pch = 20, col = "#0000ff22")

#The hemoglobin cluster is excluded based on UMI count, this is not necessarily what we want at this point. Therefore, retain UMAP_1 >5.
selected_mito <- qcdata[qcdata$mito_fraction<=0.05,] # 3536 barcodes selected_expr <- qcdata[qcdata$req_umis_per_cb>=1500,] # 2363 barcodes selected <- qcdata[qcdata$mito_fraction<=0.05 & qcdata$req_umis_per_cb >= 1500,] # 2299 barcodes
selected_blood <- qcdata[qcdata$UMAP_1 >5,] # 376 barcodes selected <- qcdata[qcdata$mito_fraction<=0.05 & (qcdata$req_umis_per_cb >= 1500 | qcdata$UMAP_1 > 5),] # 2644 barcodes
plot(log10(qcdata$mito_fraction), log10(qcdata$req_umis_per_cb), pch = 20, col = "#00000022")
points(log10(selected$mito_fraction), log10(selected$req_umis_per_cb), pch = 20, col = "#00ff0055")
write.table(selected, file ="sample1_postqc_25112019.txt", na = "0", row.names = T, sep = "\t")

