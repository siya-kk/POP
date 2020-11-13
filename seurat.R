library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(RColorBrewer)
library(DoubletFinder)


############################## Doublet removal
N1.data <- Read10X("./N1/outs/filtered_gene_bc_matrices/hg19");
colnames(x = N1.data) <- paste('N1', colnames(x = N1.data), sep = '-');
N1 <- CreateSeuratObject(counts = N1.data, min.cells = 5)
N1[["percent.mt"]] <- PercentageFeatureSet(object = N1, pattern = "MT-")
N1 <- subset(x = N1, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
N1 <- NormalizeData(object = N1, verbose = FALSE)
N1 <- FindVariableFeatures(object = N1, selection.method = "vst", nfeatures = 2000)
N1$sample <- "normal"
N1$tech <- "N1"
combined.1 <- ScaleData(object = N1, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_N1.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0253*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.21, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
N1_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.21_76=="Singlet",]), do.clean=T)



N2.data <- Read10X("./N2/outs/filtered_gene_bc_matrices/hg19");
colnames(x = N2.data) <- paste('N2', colnames(x = N2.data), sep = '-');
N2 <- CreateSeuratObject(counts = N2.data, min.cells = 5)
N2[["percent.mt"]] <- PercentageFeatureSet(object = N2, pattern = "MT-")
N2 <- subset(x = N2, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
N2 <- NormalizeData(object = N2, verbose = FALSE)
N2 <- FindVariableFeatures(object = N2, selection.method = "vst", nfeatures = 2000)
N2$sample <- "normal"
N2$tech <- "N2"
combined.1 <- ScaleData(object = N2, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_N2.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.0315*length(combined.6@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.7P@meta.data)
N2_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.06_120=="Singlet",]), do.clean=T)



N3.data <- Read10X("./N3/outs/filtered_gene_bc_matrices/hg19");
colnames(x = N3.data) <- paste('N3', colnames(x = N3.data), sep = '-');
N3 <- CreateSeuratObject(counts = N3.data, min.cells = 5)
N3[["percent.mt"]] <- PercentageFeatureSet(object = N3, pattern = "MT-")
N3 <- subset(x = N3, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
N3 <- NormalizeData(object = N3, verbose = FALSE)
N3 <- FindVariableFeatures(object = N3, selection.method = "vst", nfeatures = 2000)
N3$sample <- "normal"
N3$tech <- "N3"
combined.1 <- ScaleData(object = N3, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_N3.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0257*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.03, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
N3_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.03_73=="Singlet",]), do.clean=T)




N4.data <- Read10X("./N4/outs/filtered_gene_bc_matrices/hg19");
colnames(x = N4.data) <- paste('N4', colnames(x = N4.data), sep = '-');
N4 <- CreateSeuratObject(counts = N4.data, min.cells = 5)
N4[["percent.mt"]] <- PercentageFeatureSet(object = N4, pattern = "MT-")
N4 <- subset(x = N4, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
N4 <- NormalizeData(object = N4, verbose = FALSE)
N4 <- FindVariableFeatures(object = N4, selection.method = "vst", nfeatures = 2000)
N4$sample <- "normal"
N4$tech <- "N4"
combined.1 <- ScaleData(object = N4, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_N4.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0184*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
N4_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.06_43=="Singlet",]), do.clean=T)
dim(N4_final)





N5.data <- Read10X("./N5/outs/filtered_gene_bc_matrices/hg19");
colnames(x = N5.data) <- paste('N5', colnames(x = N5.data), sep = '-');
N5 <- CreateSeuratObject(counts = N5.data, min.cells = 5)
N5[["percent.mt"]] <- PercentageFeatureSet(object = N5, pattern = "MT-")
N5 <- subset(x = N5, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
N5 <- NormalizeData(object = N5, verbose = FALSE)
N5 <- FindVariableFeatures(object = N5, selection.method = "vst", nfeatures = 2000)
N5$sample <- "normal"
N5$tech <- "N5"
combined.1 <- ScaleData(object = N5, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_N5.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0331*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
N5_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.06_135=="Singlet",]), do.clean=T)
dim(combined.10)
dim(N5_final)





P1.data <- Read10X("./P1/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P1.data) <- paste('P1', colnames(x = P1.data), sep = '-');
P1 <- CreateSeuratObject(counts = P1.data, min.cells = 5)
P1[["percent.mt"]] <- PercentageFeatureSet(object = P1, pattern = "MT-")
P1 <- subset(x = P1, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P1 <- NormalizeData(object = P1, verbose = FALSE)
P1 <- FindVariableFeatures(object = P1, selection.method = "vst", nfeatures = 2000)
P1$sample <- "POP"
P1$tech <- "P1"
combined.1 <- ScaleData(object = P1, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P1.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0316*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.08, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P1_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.08_121=="Singlet",]), do.clean=T)




P2.data <- Read10X("./P2/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P2.data) <- paste('P2', colnames(x = P2.data), sep = '-');
P2 <- CreateSeuratObject(counts = P2.data, min.cells = 5)
P2[["percent.mt"]] <- PercentageFeatureSet(object = P2, pattern = "MT-")
P2 <- subset(x = P2, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P2 <- NormalizeData(object = P2, verbose = FALSE)
P2 <- FindVariableFeatures(object = P2, selection.method = "vst", nfeatures = 2000)
P2$sample <- "POP"
P2$tech <- "P2"
combined.1 <- ScaleData(object = P2, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P2.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0257*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P2_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.01_79=="Singlet",]), do.clean=T)




P3.data <- Read10X("./P3/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P3.data) <- paste('P3', colnames(x = P3.data), sep = '-');
P3 <- CreateSeuratObject(counts = P3.data, min.cells = 5)
P3[["percent.mt"]] <- PercentageFeatureSet(object = P3, pattern = "MT-")
P3 <- subset(x = P3, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P3 <- NormalizeData(object = P3, verbose = FALSE)
P3 <- FindVariableFeatures(object = P3, selection.method = "vst", nfeatures = 2000)
P3$sample <- "POP"
P3$tech <- "P3"
combined.1 <- ScaleData(object = P3, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P3.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0153*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.04, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P3_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.04_23=="Singlet",]), do.clean=T)





P4.data <- Read10X("./P4/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P4.data) <- paste('P4', colnames(x = P4.data), sep = '-');
P4 <- CreateSeuratObject(counts = P4.data, min.cells = 5)
P4[["percent.mt"]] <- PercentageFeatureSet(object = P4, pattern = "MT-")
P4 <- subset(x = P4, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P4 <- NormalizeData(object = P4, verbose = FALSE)
P4 <- FindVariableFeatures(object = P4, selection.method = "vst", nfeatures = 2000)
P4$sample <- "POP"
P4$tech <- "P4"
combined.1 <- ScaleData(object = P4, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P4.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0433*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.04, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P4_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.04_229=="Singlet",]), do.clean=T)



P5.data <- Read10X("./P5/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P5.data) <- paste('P5', colnames(x = P5.data), sep = '-');
P5 <- CreateSeuratObject(counts = P5.data, min.cells = 5)
P5[["percent.mt"]] <- PercentageFeatureSet(object = P5, pattern = "MT-")
P5 <- subset(x = P5, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P5 <- NormalizeData(object = P5, verbose = FALSE)
P5 <- FindVariableFeatures(object = P5, selection.method = "vst", nfeatures = 2000)
P5$sample <- "POP"
P5$tech <- "P5"
combined.1 <- ScaleData(object = P5, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P5.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0315*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.11, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P5_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.11_122=="Singlet",]), do.clean=T)




P6.data <- Read10X("./P6/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P6.data) <- paste('P6', colnames(x = P6.data), sep = '-');
P6 <- CreateSeuratObject(counts = P6.data, min.cells = 5)
P6[["percent.mt"]] <- PercentageFeatureSet(object = P6, pattern = "MT-")
P6 <- subset(x = P6, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P6 <- NormalizeData(object = P6, verbose = FALSE)
P6 <- FindVariableFeatures(object = P6, selection.method = "vst", nfeatures = 2000)
P6$sample <- "POP"
P6$tech <- "P6"
combined.1 <- ScaleData(object = P6, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P6.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0274*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.11, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P6_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.11_122=="Singlet",]), do.clean=T)




P7.data <- Read10X("./P7/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P7.data) <- paste('P7', colnames(x = P7.data), sep = '-');
P7 <- CreateSeuratObject(counts = P7.data, min.cells = 5)
P7[["percent.mt"]] <- PercentageFeatureSet(object = P7, pattern = "MT-")
P7 <- subset(x = P7, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P7 <- NormalizeData(object = P7, verbose = FALSE)
P7 <- FindVariableFeatures(object = P7, selection.method = "vst", nfeatures = 2000)
P7$sample <- "POP"
P7$tech <- "P7"
combined.1 <- ScaleData(object = P7, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P7.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0402*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P7_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.01_88=="Singlet",]), do.clean=T)





P8.data <- Read10X("./P8/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P8.data) <- paste('P8', colnames(x = P8.data), sep = '-');
P8 <- CreateSeuratObject(counts = P8.data, min.cells = 5)
P8[["percent.mt"]] <- PercentageFeatureSet(object = P8, pattern = "MT-")
P8 <- subset(x = P8, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P8 <- NormalizeData(object = P8, verbose = FALSE)
P8 <- FindVariableFeatures(object = P8, selection.method = "vst", nfeatures = 2000)
P8$sample <- "POP"
P8$tech <- "P8"
combined.1 <- ScaleData(object = P8, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P8.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0248*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.1, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P8_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.1_184=="Singlet",]), do.clean=T)




P9.data <- Read10X("./P9/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P9.data) <- paste('P9', colnames(x = P9.data), sep = '-');
P9 <- CreateSeuratObject(counts = P9.data, min.cells = 5)
P9[["percent.mt"]] <- PercentageFeatureSet(object = P9, pattern = "MT-")
P9 <- subset(x = P9, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P9 <- NormalizeData(object = P9, verbose = FALSE)
P9 <- FindVariableFeatures(object = P9, selection.method = "vst", nfeatures = 2000)
P9$sample <- "POP"
P9$tech <- "P9"
combined.1 <- ScaleData(object = P9, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P9.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0355*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.02, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P9_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.02_74=="Singlet",]), do.clean=T)




P10.data <- Read10X("./P10/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P10.data) <- paste('P10', colnames(x = P10.data), sep = '-');
P10 <- CreateSeuratObject(counts = P10.data, min.cells = 5)
P10[["percent.mt"]] <- PercentageFeatureSet(object = P10, pattern = "MT-")
P10 <- subset(x = P10, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P10 <- NormalizeData(object = P10, verbose = FALSE)
P10 <- FindVariableFeatures(object = P10, selection.method = "vst", nfeatures = 2000)
P10$sample <- "POP"
P10$tech <- "P10"
combined.1 <- ScaleData(object = P10, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P10.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0563*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.005, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P10_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.005_393=="Singlet",]), do.clean=T)





P11.data <- Read10X("./P11/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P11.data) <- paste('P11', colnames(x = P11.data), sep = '-');
P11 <- CreateSeuratObject(counts = P11.data, min.cells = 5)
P11[["percent.mt"]] <- PercentageFeatureSet(object = P11, pattern = "MT-")
P11 <- subset(x = P11, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P11 <- NormalizeData(object = P11, verbose = FALSE)
P11 <- FindVariableFeatures(object = P11, selection.method = "vst", nfeatures = 2000)
P11$sample <- "POP"
P11$tech <- "P11"
combined.1 <- ScaleData(object = P11, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P11.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0407*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P11_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.09_209=="Singlet",]), do.clean=T)




P12.data <- Read10X("./P12/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P12.data) <- paste('P12', colnames(x = P12.data), sep = '-');
P12 <- CreateSeuratObject(counts = P12.data, min.cells = 5)
P12[["percent.mt"]] <- PercentageFeatureSet(object = P12, pattern = "MT-")
P12 <- subset(x = P12, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P12 <- NormalizeData(object = P12, verbose = FALSE)
P12 <- FindVariableFeatures(object = P12, selection.method = "vst", nfeatures = 2000)
P12$sample <- "POP"
P12$tech <- "P12"
combined.1 <- ScaleData(object = P12, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P12.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0405*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.03, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P12_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.03_206=="Singlet",]), do.clean=T)




P13.data <- Read10X("./P13/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P13.data) <- paste('P13', colnames(x = P13.data), sep = '-');
P13 <- CreateSeuratObject(counts = P13.data, min.cells = 5)
P13[["percent.mt"]] <- PercentageFeatureSet(object = P13, pattern = "MT-")
P13 <- subset(x = P13, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P13 <- NormalizeData(object = P13, verbose = FALSE)
P13 <- FindVariableFeatures(object = P13, selection.method = "vst", nfeatures = 2000)
P13$sample <- "POP"
P13$tech <- "P13"
combined.1 <- ScaleData(object = P13, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P13.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0173*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P13_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.01_34=="Singlet",]), do.clean=T)





P14.data <- Read10X("./P14/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P14.data) <- paste('P14', colnames(x = P14.data), sep = '-');
P14 <- CreateSeuratObject(counts = P14.data, min.cells = 5)
P14[["percent.mt"]] <- PercentageFeatureSet(object = P14, pattern = "MT-")
P14 <- subset(x = P14, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P14 <- NormalizeData(object = P14, verbose = FALSE)
P14 <- FindVariableFeatures(object = P14, selection.method = "vst", nfeatures = 2000)
P14$sample <- "POP"
P14$tech <- "P14"
combined.1 <- ScaleData(object = P14, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P14.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0382*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P14_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.06_181=="Singlet",]), do.clean=T)





P15.data <- Read10X("./P15/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P15.data) <- paste('P15', colnames(x = P15.data), sep = '-');
P15 <- CreateSeuratObject(counts = P15.data, min.cells = 5)
P15[["percent.mt"]] <- PercentageFeatureSet(object = P15, pattern = "MT-")
P15 <- subset(x = P15, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P15 <- NormalizeData(object = P15, verbose = FALSE)
P15 <- FindVariableFeatures(object = P15, selection.method = "vst", nfeatures = 2000)
P15$sample <- "POP"
P15$tech <- "P15"
combined.1 <- ScaleData(object = P15, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P15.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0489*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.04, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P15_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.04_293=="Singlet",]), do.clean=T)




P16.data <- Read10X("./P16/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P16.data) <- paste('P16', colnames(x = P16.data), sep = '-');
P16 <- CreateSeuratObject(counts = P16.data, min.cells = 5)
P16[["percent.mt"]] <- PercentageFeatureSet(object = P16, pattern = "MT-")
P16 <- subset(x = P16, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P16 <- NormalizeData(object = P16, verbose = FALSE)
P16 <- FindVariableFeatures(object = P16, selection.method = "vst", nfeatures = 2000)
P16$sample <- "POP"
P16$tech <- "P16"
combined.1 <- ScaleData(object = P16, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P16.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0462*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.04, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P16_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.04_268=="Singlet",]), do.clean=T)





##############################Perform integration
anchors <- FindIntegrationAnchors(object.list = list(N1_final,N2_final,N3_final,N4_final,N5_final,P1_final,P2_final,P3_final,P4_final,P5_final,P6_final,P7_final,P8_final,P9_final,P10_final,P11_final,P12_final,P13_final,P14_final,P15_final,P16_final), dims = 1:40)
combined <- IntegrateData(anchorset = anchors, dims = 1:40)
DefaultAssay(object = combined) <- "integrated"



##############################Perform an integrated analysis
combined.1 <- ScaleData(object = combined, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 40, verbose = FALSE)
save(combined.2,file="RunPCA.Robj") 
pdf("DimHeatmap.pdf")
DimHeatmap(combined.2, dims = 1:20, cells = 500, balanced = TRUE)
DimHeatmap(combined.2, dims = 21:40, cells = 500, balanced = TRUE)
dev.off()
pdf("ElbowPlot.pdf")
ElbowPlot(combined.2, ndims = 40)
dev.off()
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:40)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.5 <- RunUMAP(object = combined.4, reduction = "pca", dims = 1:40)
save(combined.5,file="RunUMAP.Robj") 
pdf("DimPlot_v1.pdf")
DimPlot(combined.5, reduction = "umap", pt.size = 0.1, label=TRUE)
dev.off()



##############################Identify conserved cell type markers
all.markers <- FindAllMarkers(combined.5, logfc.threshold = 0.25)
head(all.markers)

##############################Merge & Rename
combined.6 <- SubsetData(combined.5, ident.use = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22), do.clean = TRUE)
combined.7 <- RenameIdents(combined.6, 
'0'='FIB',
'1'='FIB',
'2'='FIB',
'3'='FIB',
'4'='SMC',
'5'='EC',
'6'='SMC',
'7'='MΦ',
'8'='MEP',
'9'='MΦ',
'10'='FIB',
'11'='TC',
'12'='SMC',
'13'='MΦ',
'14'='TC',
'15'='PB',
'16'='EP',
'17'='FIB',
'18'='FIB',
'19'='LEC',
'20'='BC',
'21'='FIB',
'22'='MAST');
levels(combined.7)
pdf("DimPlot_v2.pdf",width=10,height=10)
DimPlot(combined.7, reduction = "umap", label = FALSE, pt.size = 0.8)
dev.off()
save(combined.7,file="RunUMAP_rename.Robj") 


