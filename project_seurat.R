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

# Find all markers of 8 clusters
for (i in 0:7){
  cluster.markers <- FindMarkers(pbmc, ident.1 = i, min.pct = 0.25)
  print(head(cluster.markers, n = 5))
}

# Find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
topgenes = pbmc.markers %>% group_by(cluster) %>% slice_max(n = 3, order_by = avg_log2FC)
topgenes$gene

# Gene expression analysis by cluster
FeaturePlot(pbmc, features = c(topgenes$gene[1:3]))
FeaturePlot(pbmc, features = c(topgenes$gene[4:6]))
FeaturePlot(pbmc, features = c(topgenes$gene[7:9]))
FeaturePlot(pbmc, features = c(topgenes$gene[10:12]))
FeaturePlot(pbmc, features = c(topgenes$gene[13:15]))
FeaturePlot(pbmc, features = c(topgenes$gene[16:18]))
FeaturePlot(pbmc, features = c(topgenes$gene[19:21]))
FeaturePlot(pbmc, features = c(topgenes$gene[22:24]))

# Heat maps 
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# Assigning cell type identity to clusters
new.cluster.ids <- c("Microglia", "Cholinergic neurons", "Mature oligodendrocytes", "Vascular endothelial cells", "Noradrenergic neurons", "Excitatory neurons, thalamus", 
                     "Excitatory neurons, midbrain", "Schwann cells")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# Violin plot with 3 cell types specific genes
# Cluster 0
VlnPlot(pbmc, features = c("C1qa", "Hexb", "Cst3"), slot = "counts", log = TRUE)
# Cluster 1
VlnPlot(pbmc, features = c("Neurod1", "Gm2694", "Adcy1"), slot = "counts", log = TRUE)
# Cluster 2
VlnPlot(pbmc, features = c("Cldn11", "Mog", "Plp1"), slot = "counts", log = TRUE)
# Cluster 3
VlnPlot(pbmc, features = c("Acta2", "Myl9", "Crip1"), slot = "counts", log = TRUE)
# Cluster 4
VlnPlot(pbmc, features = c("Npy", "Gal", "Prph"), slot = "counts", log = TRUE)
# Cluster 5
VlnPlot(pbmc, features = c("Nptxr", "Ntng1", "Rab3c"), slot = "counts", log = TRUE)
# Cluster 6
VlnPlot(pbmc, features = c("Nrgn", "6330403K07Rik", "Pcp4"), slot = "counts", log = TRUE)
# Cluster 7
VlnPlot(pbmc, features = c("Mpz", "Ncmap", "Pmp22"), slot = "counts", log = TRUE)


#save
saveRDS(pbmc, file = "../output/pca_mouse_10.rds")
