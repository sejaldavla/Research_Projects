# Load data
Female1 <- read.delim("/Users/sejaldavla/Desktop/Data_Portfolio/RNAseq/GSM4213596_female_rep1_dge.txt")
str(Female1)


Female2 <- read.delim("~/Documents/RNAseq data/Adult VNC/GSM4213597_female_rep2_dge.txt")
Male1 <- read.delim("~/Documents/RNAseq data/Adult VNC/GSM4213598_male_rep1_dge.txt")
Male2 <- read.delim("~/Documents/RNAseq data/Adult VNC/GSM4213599_male_rep2_dge.txt")

#Matrix is character vector, convert to numeric, assign row names as genes
row.names(Female1) <- Female1[,1]
Female1_t <- as.matrix(Female1[,-1])

# Summary counts for genes and cells
counts_per_cell <- Matrix::colSums(Female1_t, na.rm = FALSE)
counts_per_cell
mean(counts_per_cell)

counts_per_gene <- Matrix::rowSums(Female1_t, na.rm = FALSE)
counts_per_gene
mean(counts_per_gene)

genes_per_cell <- Matrix::colSums(Female1_t > 0)
genes_per_cell

cells_per_gene <- Matrix::rowSums(Female1 [,-1] > 0)
cells_per_gene
mean(cells_per_gene)

hist(log10(genes_per_cell+1), col = 'wheat', xlab = 'log10 genes per cell', ylab = 'no of cells')
hist(log10(cells_per_gene+1), col = 'wheat', xlab = 'log10 cells per gene', ylab = 'no of genes')
plot(sort(genes_per_cell), xlab = 'cell', log = 'y', main = 'genes per cell (ordered)')

# library
library(Seurat)
library(dplyr)
library(SingleCellExperiment)
library(scran)

#Create Seurat object
Female1vnc <- CreateSeuratObject(Female1_t)
Female1vnc@meta.data[1:5,1:3]
VlnPlot(Female1vnc, c("nCount_RNA","nFeature_RNA"), pt.size = 0.1)
VlnPlot(Female1vnc, "nCount_RNA", y.max = 25000, pt.size = 0.1)

#nFeature_RNA: The no of genes, nCount_RNA: the number of reads
FeatureScatter(Female1vnc, feature1 = 'nCount_RNA',feature2 = 'nFeature_RNA', xlab ='no of counts', ylab='no of genes') + NoLegend()

#Control for mitochondrial genes
mito.genes <- grep(pattern = "^mt:", x = rownames(x = Female1_t), value = TRUE)
mito.genes
percent.mito <- Matrix::colSums(Female1_t[mito.genes,]) / Matrix::colSums(Female1_t[,])
percent.mito
Female1vnc <- AddMetaData(object = Female1vnc, metadata = percent.mito, col.name = "percent.mito") # Does not work
head(Female1vnc@meta.data)

Female1vnc[["percent.mito"]] <- PercentageFeatureSet(Female1vnc, pattern = "^mt:")

# Feature plots
VlnPlot(Female1vnc, features = c("nCount_RNA","nFeature_RNA","percent.mito"), ncol = 3, pt.size = 0.1)
VlnPlot(Female1vnc, c("nCount_RNA","nFeature_RNA","percent.mito"), pt.size = 0.1)

#Normalize data
Female1vnc <- NormalizeData(Female1vnc)
#Scale data
Female1vnc <- ScaleData(Female1vnc, features = NULL, do.scale = TRUE)

Female1vnc <- FindVariableFeatures(Female1vnc, x.low.cutoff = 0.001, x.high.cutoff = Inf, y.cutoff = 0.001)
Female1vnc <- FindNeighbors(Female1vnc)
Female1vnc <- FindClusters(Female1vnc, resolution = 12)

#PCA, TSNE
Female1vnc <- RunPCA(Female1vnc, features = NULL)
Female1vnc <- RunTSNE(Female1vnc,reduction = "pca", dims = 1:50, tsne.method = "Rtsne", reduction.name = "tsne", check_duplicates = FALSE)
DimPlot(Female1vnc, reduction = 'tsne', label = FALSE) + NoLegend()
table(Idents(Female1vnc))

#UMAP
Female1vnc <- RunUMAP(Female1vnc, dims = 1:50, reduction = "pca")
DimPlot(Female1vnc, reduction = 'umap', label = FALSE) + NoLegend()

#Featureplots of important genes
FeaturePlot(Female1vnc, features = c('nSyb','per'), dims = c(1,2))
FeaturePlot(Female1vnc, features = c('Gapdh1'), dims = c(1,2))
FeaturePlot(Female1vnc, features = c('RpL11'), dims = c(1,2), pt.size = 0.5)
FeaturePlot(Female1vnc, features = c('Pdp1'), dims = c(1,2), pt.size = 0.5)
FeaturePlot(Female1vnc, features = c('vri'), dims = c(1,2), pt.size = 0.5)
FeaturePlot(Female1vnc, features = c('cry'), dims = c(1,2), pt.size = 0.5)
FeaturePlot(Female1vnc, features = c('Clk'), dims = c(1,2), pt.size = 0.5)
FeaturePlot(Female1vnc, features = c('tim'), dims = c(1,2), pt.size = 0.5)
FeaturePlot(Female1vnc, features = c('dco'), dims = c(1,2), pt.size = 0.5)

# Which clusters expressed genes of interest
VlnPlot(Female1vnc, features = 'Eaat1', pt.size = 0.1) + NoLegend()
VlnPlot(Female1vnc, features = 'Clk', pt.size = 0.1, sort = TRUE, idents = 1:10) + NoLegend()
VlnPlot(Female1vnc, features = 'tim', pt.size = 1, sort = TRUE, idents = 1:10) + NoLegend()
VlnPlot(Female1vnc, features = 'vri', pt.size = 0.1, sort = TRUE, idents = NULL) + NoLegend()
VlnPlot(Female1vnc, features = 'Pdp1', pt.size = 0.1, sort = TRUE, idents = 1:20) + NoLegend()
VlnPlot(Female1vnc, features = 'per', pt.size = 0.1, sort = TRUE, idents = 100:120) + NoLegend()
VlnPlot(Female1vnc, features = 'dco', pt.size = 0.1, sort = TRUE, idents = NULL) + NoLegend()
VlnPlot(Female1vnc, features = 'cry', pt.size = 0.1, sort = TRUE, idents = NULL) + NoLegend()
VlnPlot(Female1vnc, features = 'cwo', pt.size = 0.1, sort = TRUE, idents = NULL) + NoLegend()
VlnPlot(Female1vnc, features = 'Pdf', pt.size = 0.1, sort = TRUE, idents = 1:20) + NoLegend()
VlnPlot(Female1vnc, features = 'CCHa1', pt.size = 0.1, sort = TRUE, idents = 1:20) + NoLegend()
VlnPlot(Female1vnc, features = 'Trissin', pt.size = 0.1, sort = TRUE, idents = 1:20) + NoLegend()



#Correlation between genes/features
FeatureScatter(Female1vnc, 'nSyb', 'repo') + NoLegend()
FeatureScatter(Female1vnc, 'nFeature_RNA', 'nCount_RNA', pt.size = 1) + NoLegend() + xlab ('no of genes') + ylab ('no of reads')
FeatureScatter(Female1vnc, 'tim', 'Pdp1', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'tim', 'nSyb', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'tim', 'repo', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'Clk', 'repo', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'Pdp1', 'nSyb', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'Pdp1', 'repo', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'per', 'nSyb', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'per', 'repo', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'tim', 'Pdf', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'tim', 'VGlut', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'tim', 'Gad1', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'VGlut', 'Gad1', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'dco', 'nSyb', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'dco', 'repo', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'cry', 'nSyb', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'cry', 'repo', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'cwo', 'nSyb', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'cwo', 'repo', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'tim', 'Pdf', pt.size = 1, label = TRUE) + NoLegend()
FeatureScatter(Female1vnc, 'per', 'Pdf', pt.size = 1, label = TRUE) + NoLegend()

RidgePlot(Female1vnc, 'per') + NoLegend()
Seurat::BuildClusterTree(Female1vnc)

#cluster info
plot(Idents(Female1vnc), xlab = 'cluster', ylab = 'no of cells', col = 'cyan')

#ClusterHeatmap library
library(devtools)
library(usethis)
install_github("jokergoo/ComplexHeatmap")

#load
library(ComplexHeatmap)

#Find markers
FM <- FindAllMarkers(Female1vnc, test.use = 'negbinom', only.pos = TRUE)

#Find markers of individual clusters with different criteria
cluster1markers <- FindMarkers(Female1vnc, ident.1 = 1, min.pct = 0.25)
cluster1markers

cluster2markers <- FindMarkers(Female1vnc, ident.1 = 2, only.pos = TRUE)
cluster2markers

cluster76markers <- FindMarkers(Female1vnc, ident.1 = 76, only.pos = TRUE)
cluster76markers

cluster7markers <- FindMarkers(Female1vnc, ident.1 = 7, min.pct = 0.25)
cluster7markers

DotPlot(cluster7markers, features = row.names(cluster7markers[1:10]))

#Define top 10 based on logFC
top10 <- FM %>% group_by(cluster) %>% top_n(10, avg_logFC)
top10

DoHeatmap(Female1vnc,group.by = 'ident', disp.min = 4, disp.max = 6, slot = 'scale.data', label = TRUE)


FM <- FM %>% 
  filter(avg_logFC>0.5) %>% 
  filter(p_val_adj<0.05)

#Save as CSV
write.csv(FM,'ClusterMarkers.csv')


