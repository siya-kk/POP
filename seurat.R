library(Seurat)
library(dplyr)
library(Matrix)
library(methods)
library(RColorBrewer)
library(DoubletFinder)


##############################filter doublet
N01.data <- Read10X("./N01/outs/filtered_gene_bc_matrices/hg19");
colnames(x = N01.data) <- paste('N01', colnames(x = N01.data), sep = '-');
N01 <- CreateSeuratObject(counts = N01.data, min.cells = 5)
N01[["percent.mt"]] <- PercentageFeatureSet(object = N01, pattern = "MT-")
N01 <- subset(x = N01, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
N01 <- NormalizeData(object = N01, verbose = FALSE)
N01 <- FindVariableFeatures(object = N01, selection.method = "vst", nfeatures = 2000)
N01$sample <- "normal"
N01$tech <- "N01"
combined.1 <- ScaleData(object = N01, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_N01.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0253*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.21, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
N01_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.21_76=="Singlet",]), do.clean=T)



N02.data <- Read10X("./N0/outs/filtered_gene_bc_matrices/hg19");
colnames(x = N02.data) <- paste('N02', colnames(x = N02.data), sep = '-');
N02 <- CreateSeuratObject(counts = N02.data, min.cells = 5)
N02[["percent.mt"]] <- PercentageFeatureSet(object = N02, pattern = "MT-")
N02 <- subset(x = N02, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
N02 <- NormalizeData(object = N02, verbose = FALSE)
N02 <- FindVariableFeatures(object = N02, selection.method = "vst", nfeatures = 2000)
N02$sample <- "normal"
N02$tech <- "N02"
combined.1 <- ScaleData(object = N02, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_N02.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.0315*length(combined.6@active.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.7P@meta.data)
N02_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.06_120=="Singlet",]), do.clean=T)



N03.data <- Read10X("./N03/outs/filtered_gene_bc_matrices/hg19");
colnames(x = N03.data) <- paste('N03', colnames(x = N03.data), sep = '-');
N03 <- CreateSeuratObject(counts = N03.data, min.cells = 5)
N03[["percent.mt"]] <- PercentageFeatureSet(object = N03, pattern = "MT-")
N03 <- subset(x = N03, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
N03 <- NormalizeData(object = N03, verbose = FALSE)
N03 <- FindVariableFeatures(object = N03, selection.method = "vst", nfeatures = 2000)
N03$sample <- "normal"
N03$tech <- "N03"
combined.1 <- ScaleData(object = N03, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_N03.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0257*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.03, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
N03_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.03_73=="Singlet",]), do.clean=T)




N04.data <- Read10X("./N04/outs/filtered_gene_bc_matrices/hg19");
colnames(x = N04.data) <- paste('N04', colnames(x = N04.data), sep = '-');
N04 <- CreateSeuratObject(counts = N04.data, min.cells = 5)
N04[["percent.mt"]] <- PercentageFeatureSet(object = N04, pattern = "MT-")
N04 <- subset(x = N04, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
N04 <- NormalizeData(object = N04, verbose = FALSE)
N04 <- FindVariableFeatures(object = N04, selection.method = "vst", nfeatures = 2000)
N04$sample <- "normal"
N04$tech <- "N04"
combined.1 <- ScaleData(object = N04, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_N04.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0184*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
N04_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.06_43=="Singlet",]), do.clean=T)
dim(N04_final)





N05.data <- Read10X("./N05/outs/filtered_gene_bc_matrices/hg19");
colnames(x = N05.data) <- paste('N05', colnames(x = N05.data), sep = '-');
N05 <- CreateSeuratObject(counts = N05.data, min.cells = 5)
N05[["percent.mt"]] <- PercentageFeatureSet(object = N05, pattern = "MT-")
N05 <- subset(x = N05, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
N05 <- NormalizeData(object = N05, verbose = FALSE)
N05 <- FindVariableFeatures(object = N05, selection.method = "vst", nfeatures = 2000)
N05$sample <- "normal"
N05$tech <- "N05"
combined.1 <- ScaleData(object = N05, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_N05.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0331*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.06, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
N05_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.06_135=="Singlet",]), do.clean=T)
dim(combined.10)
dim(N05_final)





P01.data <- Read10X("./P01/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P01.data) <- paste('P01', colnames(x = P01.data), sep = '-');
P01 <- CreateSeuratObject(counts = P01.data, min.cells = 5)
P01[["percent.mt"]] <- PercentageFeatureSet(object = P01, pattern = "MT-")
P01 <- subset(x = P01, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P01 <- NormalizeData(object = P01, verbose = FALSE)
P01 <- FindVariableFeatures(object = P01, selection.method = "vst", nfeatures = 2000)
P01$sample <- "POP"
P01$tech <- "P01"
combined.1 <- ScaleData(object = P01, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P01.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0316*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.08, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P01_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.08_121=="Singlet",]), do.clean=T)




P02.data <- Read10X("./P02/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P02.data) <- paste('P02', colnames(x = P02.data), sep = '-');
P02 <- CreateSeuratObject(counts = P02.data, min.cells = 5)
P02[["percent.mt"]] <- PercentageFeatureSet(object = P02, pattern = "MT-")
P02 <- subset(x = P02, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P02 <- NormalizeData(object = P02, verbose = FALSE)
P02 <- FindVariableFeatures(object = P02, selection.method = "vst", nfeatures = 2000)
P02$sample <- "POP"
P02$tech <- "P02"
combined.1 <- ScaleData(object = P02, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P02.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0257*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P02_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.01_79=="Singlet",]), do.clean=T)




P03.data <- Read10X("./P03/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P03.data) <- paste('P03', colnames(x = P03.data), sep = '-');
P03 <- CreateSeuratObject(counts = P03.data, min.cells = 5)
P03[["percent.mt"]] <- PercentageFeatureSet(object = P03, pattern = "MT-")
P03 <- subset(x = P03, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P03 <- NormalizeData(object = P03, verbose = FALSE)
P03 <- FindVariableFeatures(object = P03, selection.method = "vst", nfeatures = 2000)
P03$sample <- "POP"
P03$tech <- "P03"
combined.1 <- ScaleData(object = P03, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P03.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0153*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.04, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P03_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.04_23=="Singlet",]), do.clean=T)





P04.data <- Read10X("./P04/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P04.data) <- paste('P04', colnames(x = P04.data), sep = '-');
P04 <- CreateSeuratObject(counts = P04.data, min.cells = 5)
P04[["percent.mt"]] <- PercentageFeatureSet(object = P04, pattern = "MT-")
P04 <- subset(x = P04, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P04 <- NormalizeData(object = P04, verbose = FALSE)
P04 <- FindVariableFeatures(object = P04, selection.method = "vst", nfeatures = 2000)
P04$sample <- "POP"
P04$tech <- "P04"
combined.1 <- ScaleData(object = P04, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P04.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0433*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.04, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P04_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.04_229=="Singlet",]), do.clean=T)



P05.data <- Read10X("./P05/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P05.data) <- paste('P05', colnames(x = P05.data), sep = '-');
P05 <- CreateSeuratObject(counts = P05.data, min.cells = 5)
P05[["percent.mt"]] <- PercentageFeatureSet(object = P05, pattern = "MT-")
P05 <- subset(x = P05, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P05 <- NormalizeData(object = P05, verbose = FALSE)
P05 <- FindVariableFeatures(object = P05, selection.method = "vst", nfeatures = 2000)
P05$sample <- "POP"
P05$tech <- "P05"
combined.1 <- ScaleData(object = P05, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P05.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0315*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.11, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P05_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.11_122=="Singlet",]), do.clean=T)




P06.data <- Read10X("./P06/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P06.data) <- paste('P06', colnames(x = P06.data), sep = '-');
P06 <- CreateSeuratObject(counts = P06.data, min.cells = 5)
P06[["percent.mt"]] <- PercentageFeatureSet(object = P06, pattern = "MT-")
P06 <- subset(x = P06, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P06 <- NormalizeData(object = P06, verbose = FALSE)
P06 <- FindVariableFeatures(object = P06, selection.method = "vst", nfeatures = 2000)
P06$sample <- "POP"
P06$tech <- "P06"
combined.1 <- ScaleData(object = P06, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P06.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0274*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.11, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P06_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.11_122=="Singlet",]), do.clean=T)




P07.data <- Read10X("./P07/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P07.data) <- paste('P07', colnames(x = P07.data), sep = '-');
P07 <- CreateSeuratObject(counts = P07.data, min.cells = 5)
P07[["percent.mt"]] <- PercentageFeatureSet(object = P07, pattern = "MT-")
P07 <- subset(x = P07, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P07 <- NormalizeData(object = P07, verbose = FALSE)
P07 <- FindVariableFeatures(object = P07, selection.method = "vst", nfeatures = 2000)
P07$sample <- "POP"
P07$tech <- "P07"
combined.1 <- ScaleData(object = P07, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P07.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0402*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P07_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.01_88=="Singlet",]), do.clean=T)





P08.data <- Read10X("./P08/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P08.data) <- paste('P08', colnames(x = P08.data), sep = '-');
P08 <- CreateSeuratObject(counts = P08.data, min.cells = 5)
P08[["percent.mt"]] <- PercentageFeatureSet(object = P08, pattern = "MT-")
P08 <- subset(x = P08, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P08 <- NormalizeData(object = P08, verbose = FALSE)
P08 <- FindVariableFeatures(object = P08, selection.method = "vst", nfeatures = 2000)
P08$sample <- "POP"
P08$tech <- "P08"
combined.1 <- ScaleData(object = P08, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P08.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0248*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.1, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P08_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.1_184=="Singlet",]), do.clean=T)




P09.data <- Read10X("./P09/outs/filtered_gene_bc_matrices/hg19")
colnames(x = P09.data) <- paste('P09', colnames(x = P09.data), sep = '-');
P09 <- CreateSeuratObject(counts = P09.data, min.cells = 5)
P09[["percent.mt"]] <- PercentageFeatureSet(object = P09, pattern = "MT-")
P09 <- subset(x = P09, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
P09 <- NormalizeData(object = P09, verbose = FALSE)
P09 <- FindVariableFeatures(object = P09, selection.method = "vst", nfeatures = 2000)
P09$sample <- "POP"
P09$tech <- "P09"
combined.1 <- ScaleData(object = P09, verbose = FALSE)
combined.2 <- RunPCA(object = combined.1, npcs = 10, verbose = FALSE)
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:10)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.6 <- RunTSNE(object = combined.4, reduction = "pca", dims = 1:10)
combined.7 <- paramSweep_v3(combined.6, PCs = 1:10,sct = FALSE)
combined.8 <- summarizeSweep(combined.7, GT = FALSE)
combined.9 <- find.pK(combined.8)
pdf("pK_P09.pdf")
plot(combined.9$pK, combined.9$BCmetric, type='b')
dev.off()
annotations <- combined.6@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.0355*length(combined.6@active.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
combined.10 <- doubletFinder_v3(combined.6, PCs = 1:10, pN = 0.25, pK = 0.02, nExp = nExp_poi, reuse.pANN = FALSE)
head(combined.10@meta.data)
P09_final <- SubsetData(combined.10, cells= rownames(combined.10@meta.data[combined.10@meta.data$DF.classifications_0.25_0.02_74=="Singlet",]), do.clean=T)




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





##############################PCA
anchors <- FindIntegrationAnchors(object.list = list(N01_final,N02_final,N03_final,N04_final,N05_final,P01_final,P02_final,P03_final,P04_final,P05_final,P06_final,P07_final,P08_final,P09_final,P10_final,P11_final,P12_final,P13_final,P14_final,P15_final,P16_final), dims = 1:40)
combined <- IntegrateData(anchorset = anchors, dims = 1:40)
DefaultAssay(object = combined) <- "integrated"
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





##############################UMAP
combined.3 <- FindNeighbors(object = combined.2, reduction = "pca", dims = 1:40)
combined.4 <- FindClusters(combined.3, resolution = 0.5)
combined.5 <- RunUMAP(object = combined.4, reduction = "pca", dims = 1:40)
save(combined.5,file="RunUMAP.Robj") 
pdf("DimPlot.pdf")
DimPlot(combined.6, reduction = "umap", pt.size = 0.1, label=TRUE)
dev.off()






