library(dplyr)
library(Seurat)
library(patchwork)

library(readr)

matrix.count <-read.csv('../data/MOUSE_BRAIN_DATASET_2_COUNTS.tsv', header=TRUE, sep='\t', row.names = 'X' ) # load count matrix 
head(matrix.count[1:5]) #verify the look of the matrix

data.matrix <-CreateSeuratObject(counts = matrix.count, project = "mouse_2", min.cells = 3, min.features = 200)

head(data.matrix)
data.matrix[["percent.mt"]] <- PercentageFeatureSet(data.matrix, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(data.matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

head(data.matrix@meta.data, 5) 

#feature scatter
plot1 <- FeatureScatter(data.matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data.matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Normalizing the data 
pbmc <- NormalizeData(data.matrix)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scaling data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#peforming linear dim reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#more plots
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

#clustering the cells 
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)

#run non-linear dim reduction
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

#save
saveRDS(pbmc, file = "../output/pca_mouse_10.rds")
