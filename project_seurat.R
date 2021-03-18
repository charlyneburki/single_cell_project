library(dplyr)
library(Seurat)
library(patchwork)
library(readr)
library(ggplot2)
library(ggpubr)

# Load count matrix
#add filepath to data location below
data_location <- '../data/MOUSE_BRAIN_DATASET_2_COUNTS.tsv' 
matrix.count <-read.csv(data_location, header=TRUE, sep='\t', row.names = 'X')

# Initialize the Seurat object
data.matrix <-CreateSeuratObject(counts = matrix.count, project = "mouse_2", min.cells = 3, min.features = 200)

# Add percent mt
data.matrix[["percent.mt"]] <- PercentageFeatureSet(data.matrix, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(data.matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Feature scatter
plot1 <- FeatureScatter(data.matrix, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data.matrix, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggarrange(plot1,plot2, nrow=2)

# Select the data and remove unwanted cells 
mouse2 <- subset(data.matrix, subset = nFeature_RNA > 100 & nFeature_RNA < 4500 & percent.mt < 5)

#Normalize the data 
mouse2 <- NormalizeData(data.matrix)

# Find variable features
mouse2 <- FindVariableFeatures(mouse2, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mouse2), 10)

# Plot variable features with labels
plot1 <- VariableFeaturePlot(mouse2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggarrange(plot2,nrow = 1)

# Scale the data
all.genes <- rownames(mouse2)
mouse2 <- ScaleData(mouse2, features = all.genes)

# Perform linear dim reduction
mouse2 <- RunPCA(mouse2, features = VariableFeatures(object = mouse2))

# Visualize PCA results in different ways
VizDimLoadings(mouse2, dims = 1:2, reduction = "pca")

DimPlot(mouse2, reduction = "pca")

DimHeatmap(mouse2, dims = 1:15, cells = 500, balanced = TRUE)

# More plots
mouse2 <- JackStraw(mouse2, num.replicate = 100)
mouse2 <- ScoreJackStraw(mouse2, dims = 1:20)
plot1 <-JackStrawPlot(mouse2, dims = 1:15)
plot2 <- ElbowPlot(mouse2)
ggarrange(plot1, plot2, nrow = 2)

# CLUSTER THE CELLS WITH A RESOLUTION OF 0.5
mouse2_1 <- FindNeighbors(mouse2, dims = 1:10)
mouse2_1 <- FindClusters(mouse2_1, resolution = 0.5)

# Run non-linear dimensional reduction
mouse2_1 <- RunUMAP(mouse2_1, dims = 1:10)
DimPlot(mouse2_1, reduction = "umap")

# FIND MARKER GENES WITH THE DEFAULT WILCOXON RANK TEST
for (i in 0:7){
  cluster.markers_1 <- FindMarkers(mouse2_1, ident.1 = i, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25)
  print('New cluster')
  print(head(cluster.markers_1, n = 5))
}

# Find marker genes for each cluster compared to all remaining cells
# With the default Wilcoxon Rank Sum test
# Report only the positive ones
mouse2.markers_1 <- FindAllMarkers(mouse2_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

topgenes <- mouse2.markers_1 %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)
top10_1 <- mouse2.markers_1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top15 <- mouse2.markers_1 %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)

# Heat map of the top 10 genes
DoHeatmap(mouse2_1, features = top10_1$gene) + NoLegend()

# Visualize top 5 marker genes expression in two different ways
VlnPlot(mouse2_1, features = c("C1qa", "Ctss", "P2ry12", "C1qb", "C1qc"))
VlnPlot(mouse2_1, features = c("Neurod1", "Gm2694", "Cbln3", "Gabra6", "Cbln1"))

VlnPlot(mouse2_1, features = c("Cldn11", "Mog", "Apod", "Mag", "Ugt8a"))
VlnPlot(mouse2_1, features = c("Acta2", "Tagln", "Myh11", "Tpm2", "Igfbp7"))
VlnPlot(mouse2_1, features = c("Npy", "Gal", "Slc6a2", "Ntrk1", "Tspan8"))
VlnPlot(mouse2_1, features = c("Slc17a6", "Necab1", "Tcf7l2", "Slc1a2", "Cit"))
VlnPlot(mouse2_1, features = c("Meis2", "Nrgn", "Ly6h", "Tcf7l2", "Hap1"))
VlnPlot(mouse2_1, features = c("Mpz", "Prx", "Cldn19", "Dhh", "Egfl8"))

FeaturePlot(mouse2_1, features = c(topgenes$gene[0:5]))
FeaturePlot(mouse2_1, features = c(topgenes$gene[15:20]))
FeaturePlot(mouse2_1, features = c(topgenes$gene[30:35]))
FeaturePlot(mouse2_1, features = c(topgenes$gene[45:50]))
FeaturePlot(mouse2_1, features = c(topgenes$gene[60:65]))
FeaturePlot(mouse2_1, features = c(topgenes$gene[75:80]))
FeaturePlot(mouse2_1, features = c(topgenes$gene[90:95]))
FeaturePlot(mouse2_1, features = c(topgenes$gene[105:110]))

# FIND MARKER GENES WITH THE STUDENT T TEST --> leads to same results
for (i in 0:7){
  cluster.markers_2 <- FindMarkers(mouse2_1, ident.1 = i, min.pct = 0.1, test.use = "t", only.pos = TRUE)
  print('New cluster')
  print(head(cluster.markers_2, n = 5))
}

# Assigning cell type identity to clusters
new.cluster.ids <- c("Microglia", "Granule neurons, cerebellum", 
                     "Mature oligodendrocytes", "Vascular smooth muscle cells, arterial",
                     "Noradrenergic neurons", "Excitatory neurons,thalamus/hypothalamus",
                     "Excitatory neurons, midbrain", "Schwann cells")

names(new.cluster.ids) <- levels(mouse2_1)
mouse2_1 <- RenameIdents(mouse2_1, new.cluster.ids)
DimPlot(mouse2_1, reduction = "umap", label = TRUE, pt.size = 0.5, repel= TRUE) + NoLegend()

# CLUSTER THE CELLS WITH A RESOLUTION OF 1.2
mouse2_2 <- FindNeighbors(mouse2, dims = 1:10)
mouse2_2 <- FindClusters(mouse2_2, resolution = 1.2)

# Run non-linear dimensional reduction
mouse2_2 <- RunUMAP(mouse2_2, dims = 1:10)
DimPlot(mouse2_2, reduction = "umap")

# FIND MARKER GENES WITH THE DEFAULT WILCOXON RANK TEST
for (i in 0:13){
  cluster.markers_3 <- FindMarkers(mouse2_2, ident.1 = i, min.pct = 0.25)
  print('New cluster')
  print(head(cluster.markers_3, n = 5))
}

# Find marker genes for each cluster compared to all remaining cells
# With the default Wilcoxon Rank Sum test
# Report only the positive ones
mouse2.markers_2 <- FindAllMarkers(mouse2_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10_2 <- mouse2.markers_2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

# Heat map of the top 10 genes
DoHeatmap(mouse2_2, features = top10_2$gene) + NoLegend()


##creating the boxplots for the dataset infos
## Info dataset
# Number of cells in the initial dataset
length(data.matrix$orig.ident)

# Number of cells after quality filtering
mouse2_1 <- subset(mouse2_1, subset = nFeature_RNA > 100 & nFeature_RNA < 4500 & percent.mt < 5)
length(mouse2_1$orig.ident)

# Number of genes in the dataset
data.matrix <-CreateSeuratObject(counts = matrix.count, project = "mouse_2", min.cells = 3, min.features = 200)
# Read information in...
data.matrix

# Mean number of genes expressed per cell
mean(data.matrix@meta.data$nFeature_RNA)

# To get how many cells per cluster
summary(mouse2_1$seurat_clusters)

# To get how many genes are expressed in each cluster in mean
mean_nb_genes_cluster <- c()
for (i in 0:7){
  mean_nb_genes_cluster <- c(mean_nb_genes_cluster, round(mean(mouse2_1$nFeature_RNA[mouse2_1$seurat_clusters==i])))
  print(round(mean(mouse2_1$nFeature_RNA[mouse2_1$seurat_clusters==i])))
}


# Create barchart with plotly
library(plotly)
fig <- plot_ly(x = c("0","1","2","3","4","5","6","7"), y = mean_nb_genes_cluster, type = "bar", text = mean_nb_genes_cluster, textposition = 'auto', name = "Mean number of genes expressed per cluster")                                 
fig %>% layout(title = "Genes expressed per cluster",
               xaxis = list(title = "Cluster"),
               yaxis = list(title = "Number of genes"))

# Boxplot
library(ggplot2)
x = mouse2_1$seurat_clusters
y =mouse2_1$nFeature_RNA
df <- tibble(x,y)
p = ggplot(df, aes(x = x, y = mouse2_1$nFeature_RNA)) + geom_jitter(shape=16, position=position_jitter(0.3), alpha=0.4) +
  geom_boxplot(outlier.shape = NA, varwidth=TRUE, color="blue")+
  scale_fill_brewer(palette="RdBu") + 
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="red")+
  stat_summary(fun=mean, geom="text", size=3, color="red", show.legend = FALSE, 
               vjust=-0.2,hjust=1.5, aes( label=round(..y.., digits=1)))+
  theme_minimal() + labs(title="Distribution of expressed genes number per cluster",x="Cluster", y = "Expressed genes")
p

