###################################
# Single-Nucleus RNA-seq Analysis #
# Pancreatic Organoids            #
# Danny Conrad                    #
# 2/15/21                         #
###################################

# Load Packages
library(Seurat)
library(ggplot2)
library(stringr)
library(caret)
library(pheatmap)
library(gridExtra)

# Set working directory
wd <- "/Volumes/DannyData/Ling_MULTIseq_1"
setwd(wd)

# Load filtered matrices from Cell Ranger output
L1.filt.data <- Read10X('outs_L1/filtered_feature_bc_matrix')
L2.filt.data <- Read10X('outs_L2/filtered_feature_bc_matrix')

# Remove uninformative genes (a.k.a. not expressed or detected in almost any nuclei)
L1.gene.counts <- Matrix::rowSums(L1.filt.data)
L1.genes <- rownames(L1.filt.data)[which(L1.gene.counts >= 3)] # minimum 3 UMIs in whole nuclei pool
L1.filt.data <- L1.filt.data[L1.genes, ]

L2.gene.counts <- Matrix::rowSums(L2.filt.data)
L2.genes <- rownames(L2.filt.data)[which(L2.gene.counts >= 3)] # minimum 3 UMIs in whole nuclei pool
L2.filt.data <- L2.filt.data[L2.genes, ]

# Remove "-1" suffix from end of cellnames
str_sub(colnames(L1.filt.data), -2, -1) <- ""
str_sub(colnames(L2.filt.data), -2, -1) <- ""

# Integrate separate lanes into one dataset
obj_list <- list(L1 = CreateSeuratObject(L1.filt.data),
                 L2 = CreateSeuratObject(L2.filt.data))
obj_list[[1]]@meta.data$Lane <- "1"
obj_list[[2]]@meta.data$Lane <- "2"
for (i in 1:length(obj_list)) {
  obj_list[[i]] <- SCTransform(obj_list[[i]], verbose = FALSE)
}
features <- SelectIntegrationFeatures(obj_list, nfeatures = 3000)
options(future.globals.maxSize = 4000 * 1024^2)
obj_list <- PrepSCTIntegration(obj_list, features, verbose = FALSE, assay = 'SCT')
anchors <- FindIntegrationAnchors(obj_list, normalization.method = "SCT", 
                                  anchor.features = features, verbose = FALSE)

# New integrated Seurat object is called "org"
org <- IntegrateData(anchorset = anchors, normalization.method = "SCT", 
                     verbose = FALSE)
org <- SCTransform(org)
org <- RunPCA(org)
org <- RunUMAP(org, dims = 1:32)
org <- FindNeighbors(org, dims = 1:32)
org <- FindClusters(org, resolution = 0.2)

# Calculate percent mitochondrial UMIs per nucleus (less informative than for whole cells because we intentionally lyse cells)
mito.genes <- grep("^MT-", rownames(org@assays$RNA@counts), value = TRUE)
percent.mito <- Matrix::colSums(org@assays$RNA@counts[mito.genes, ]) / Matrix::colSums(org@assays$RNA@counts)
org@meta.data[,"PercentMito"] <- percent.mito

# Add MULTI-seq sample classification information to Seurat (switch to "MULTIseqSampleAssignment.R" to generate these)
final.calls_1 <- readRDS("multi_barcodes/final.calls_1.rds")
final.calls_2 <- readRDS("multi_barcodes/final.calls_2.rds")
names(final.calls_1) <- paste(names(final.calls_1), "_1", sep = "") 
names(final.calls_2) <- paste(names(final.calls_2), "_2", sep = "") 

final.calls <- c(final.calls_1, final.calls_2)

org@meta.data[,"MULTI"] <- final.calls[rownames(org@meta.data)]

# Visualize MULTIseq classifications
DimPlot(org, group.by = "MULTI")

# Remove Doublets, Negatives, & Mitochondrial Genes
org_sub <- subset(org,
                  cells=rownames(org@meta.data)[which(!org@meta.data$MULTI %in% c("Doublet", "Negative"))],
                  features = rownames(org)[which(!rownames(org) %in% mito.genes)])

# Add Sample Metadata
md <- data.frame(row.names = names(bar.ref),
                 Sample = factor(c(rep("Progenitor", 2),
                                   rep("Ductal Organoid - Day 8", 2),
                                   rep("Ductal Organoid - Day 16", 2),
                                   rep("Acinar Organoid - Day 8", 3),
                                   rep("Acinar Organoid - Day 16", 2)),
                                 levels = c("Progenitor", "Ductal Organoid - Day 8", "Ductal Organoid - Day 16", "Acinar Organoid - Day 8", "Acinar Organoid - Day 16"),
                                 ordered = T),
                 Type = factor(c(rep("Progenitor", 2),
                                 rep("Ductal Organoid", 4),
                                 rep("Acinar Organoid", 5)),
                               levels = c("Progenitor", "Ductal Organoid", "Acinar Organoid"),
                               ordered = T),
                 Timepoint = factor(c(rep("Day 0", 2),
                                      rep("Day 8", 2),
                                      rep("Day 16", 2),
                                      rep("Day 8", 3),
                                      rep("Day 16", 2)),
                                    levels = c("Day 0", "Day 8", "Day 16"),
                                    ordered = T),
                 Replicate = factor(c(1,2,1,3,1,2,1,2,3,1,2),
                                    levels = 1:3,
                                    ordered = T)
)

for (i in 1:nrow(md)) {
  for (j in 1:ncol(md)) {
    org_sub@meta.data[which(org_sub@meta.data$MULTI == rownames(md)[i]),colnames(md)[j]] <- md[i,j]
  }
}

# Re-factor once in Seurat object
org_sub@meta.data$MULTI <- factor(org_sub@meta.data$MULTI,
                               levels = c("Progenitor_1", "Progenitor_2", "DuctOrg8_1", "DuctOrg8_3", "DuctOrg16_1", "DuctOrg16_2", "AcinOrg8_1", "AcinOrg8_2", "AcinOrg8_3", "AcinOrg16_1", "AcinOrg16_2"))
org_sub@meta.data$Timepoint <- factor(org_sub@meta.data$Timepoint,
                                   levels = c("Day 0", "Day 8", "Day 16"),
                                   ordered = T)
org_sub@meta.data$Sample <- factor(org_sub@meta.data$Sample,
                                levels = c("Progenitor", "Ductal Organoid - Day 8", "Ductal Organoid - Day 16", "Acinar Organoid - Day 8", "Acinar Organoid - Day 16"))
org_sub@meta.data$Replicate <- factor(org_sub@meta.data$Replicate,
                                   levels = 1:3,
                                   ordered = T)

# Remove Day 16 Organoids (poor recovery for acinar population)
org_sub <- subset(org_sub, 
                  cells=rownames(org_sub@meta.data)[which(!org_sub@meta.data$Timepoint %in% c("Day 16"))])

org_sub <- SCTransform(org_sub)
org_sub <- RunPCA(org_sub)
org_sub <- RunUMAP(org_sub, dims = 1:32)
org_sub <- FindNeighbors(org_sub, dims = 1:32)
org_sub <- FindClusters(org_sub, resolution = 0.2)

org_sub@meta.data$Cluster[which(org_sub@meta.data$seurat_clusters %in% c(1,2))] <- "Acinar Organoid"
org_sub@meta.data$Cluster[which(org_sub@meta.data$seurat_clusters %in% c(3))] <- "Ductal Organoid"
org_sub@meta.data$Cluster[which(org_sub@meta.data$seurat_clusters %in% c(0,4))] <- "Progenitor"

org_sub@meta.data$Cluster <- factor(org_sub@meta.data$Cluster,
                                    levels = c("Progenitor", "Ductal Organoid", "Acinar Organoid"),
                                    ordered = T)

# Load Single-Nucleus RNA-seq of neonatal pancreas from Tosti et al, 2020
neo <- readRDS("tosti_2020/neonatal_pancreas_2020.rds")

neo_sub <- subset(neo, 
                  cells = rownames(neo@meta.data)[which(neo@meta.data$Cluster %in% c("Acinar-s", "Acinar-i", "Ductal"))])

neo_sub <- SCTransform(neo_sub)
neo_sub <- RunPCA(neo_sub)
neo_sub <- RunUMAP(neo_sub, dims = 1:32)
neo_sub <- FindNeighbors(neo_sub, dims = 1:32)
neo_sub <- FindClusters(neo_sub, resolution = 0.2)

neo_sub@meta.data$Type[which(neo_sub$seurat_clusters %in% c(1,2,3))] <- "Acinar"
neo_sub@meta.data$Type[which(neo_sub$seurat_clusters %in% c(0))] <- "Ductal"

neo_sub <- subset(neo_sub,
                  cells = rownames(neo_sub@meta.data)[which(neo_sub@meta.data$seurat_clusters != 4)])

# Marker gene calculation
Idents(org_sub) <- org_sub@meta.data$Cluster
Idents(neo_sub) <- neo_sub@meta.data$Type

markers_org <- FindAllMarkers(org_sub, logfc.threshold = 0)
markers_neo <- FindAllMarkers(neo_sub, logfc.threshold = 0)

avg_log2FC <- 0
p_val_adj <- 0.05

markers_org <- markers_org[markers_org$p_val_adj <= p_val_adj,]
markers_neo <- markers_neo[markers_neo$p_val_adj <= p_val_adj & markers_neo$avg_log2FC > avg_log2FC,]

markers_org <- markers_org[order(markers_org$avg_log2FC, decreasing = T),]
markers_neo <- markers_neo[order(markers_neo$avg_log2FC, decreasing = T),]

# Heatmap of Markers
pdf("Final/MarkerHeatmap.pdf", width = 6, height = 22, onefile = T)
heatmap.features <- markers_org[order(markers_org$cluster),]
heatmap.features <- unique(heatmap.features$gene[heatmap.features$avg_log2FC > log(1.2)])
DoHeatmap(subset(org_sub, downsample = min(table(Idents(org_sub)))), features = heatmap.features, raster = F, disp.max = 1.5, disp.min = -1.5)
dev.off()

# Venn Diagram to show overlap of markers between organoids & progenitors
markers.pairwise.nominpct <- list(AcinarOrganoid1 = FindMarkers(org_sub, ident.1 = "Acinar Organoid", ident.2 = "Ductal Organoid", logfc.threshold = 0, min.pct = 0),
                                  AcinarOrganoid2 = FindMarkers(org_sub, ident.1 = "Acinar Organoid", ident.2 = "Progenitor", logfc.threshold = 0, min.pct = 0),
                                  DuctalOrganoid1 = FindMarkers(org_sub, ident.1 = "Ductal Organoid", ident.2 = "Acinar Organoid", logfc.threshold = 0, min.pct = 0),
                                  DuctalOrganoid2 = FindMarkers(org_sub, ident.1 = "Ductal Organoid", ident.2 = "Progenitor", logfc.threshold = 0, min.pct = 0),
                                  Progenitor1 = FindMarkers(org_sub, ident.1 = "Progenitor", ident.2 = "Acinar Organoid", logfc.threshold = 0, min.pct = 0),
                                  Progenitor2 = FindMarkers(org_sub, ident.1 = "Progenitor", ident.2 = "Ductal Organoid", logfc.threshold = 0, min.pct = 0))

markers.pairwise.nominpct <- lapply(markers.pairwise.nominpct, function(x) {
  x <- x[x$p_val_adj < 0.05 & x$avg_log2FC > 0,]
})

venn.pairwise.nominpct <- list(AcinarOrganoid = unique(c(rownames(markers.pairwise.nominpct[[1]]), rownames(markers.pairwise.nominpct[[2]]))),
                               DuctalOrganoid = unique(c(rownames(markers.pairwise.nominpct[[3]]), rownames(markers.pairwise.nominpct[[4]]))),
                               Progenitor = unique(c(rownames(markers.pairwise.nominpct[[5]]), rownames(markers.pairwise.nominpct[[6]]))))
venn.diagram(venn.pairwise.nominpct, "Final/Venn_Pairwise_NoMinPct.png", width = 600, height= 600, resolution = 100)


markers.nominpct <- FindAllMarkers(org_sub, logfc.threshold = 0, min.pct = 0, only.pos = T)
markers.nominpct <- markers.nominpct[markers.nominpct$p_val_adj < 0.05,]
write.csv(markers.nominpct, "Final/Markers_NoMinPct.csv")

# Generate table of average expression per cluster for canonical markers
av.exp <- AverageExpression(org_sub, features = c("PDX1", "SOX9", "HNF1B", "HNF1A", "RBPJL", "RBPJ", "CPA2", "CEL", "PNLIP", "CTRB1", "CTRC"))$RNA
write.csv(av.exp, "Final/AverageExpression_CanonicalMarkers.csv")


# Pearson correlation of gene expression between datasets/cell types

# Use marker gene overlap by dataset but downsample equal number of genes per cell type
# Because of random sampling, perform this 100 times and take average correlation

marker.list <- list(FindAllMarkers(org_sub, logfc.threshold = 0),
                    FindAllMarkers(neo_sub, logfc.threshold = 0))

n <- 100
res <- matrix(0L, 5, 5)
ngene <- 0
min <- 1

for (i in 1:n) {
  g.list <- marker.list
  g.list <- lapply(g.list, function(x) {
    df <- x
    df$cluster <- factor(df$cluster)
    x <- df})
  g.list <- lapply(g.list, function(x) {downSample(x, x$cluster)})
  
  # use significant markers from our dataset as genelist for comparison
  g.list[[3]] <- g.list[[1]][g.list[[1]]$gene %in% rownames(neo_sub) & g.list[[1]]$p_val_adj <= 0.05,]

  g.list <- unique(g.list[[3]]$gene)
  
  expr <- cbind(AverageExpression(org_sub, features = g.list, assays = "SCT")[[1]],
                AverageExpression(neo_sub, features = g.list, assays = "SCT")[[1]])
  
  # MALAT1 is highly represented in single-nucleus datasets and very strongly affects correlations, so we remove this one gene
  expr <- expr[grep("MALAT1", rownames(expr), invert = T),]
  
  ngene <- ngene + length(g.list)
  
  res <- res + cor(expr)
  min <- min(min, min(res))
}

res <- res/n
ngene <- ngene/n


# Show just comparisons between organoids and ductal/acinar cells
col <- colorRampPalette(c("#FFFFFF", "#D1E5F0", "#92C5DE",
                          "#4393C3", "#2166AC", "#053061")) # white -> blue

res1 <- res[4, c(-4,-5)]
res2 <- res[5, c(-5,-4)]

min <- 0.1
max <- 0.4

plot_list <- list()

p1 <- pheatmap(t(res1), display_numbers=TRUE, number_color="black", border_color=NA, cluster_cols=FALSE, cluster_rows=FALSE, 
               main="Acinar Cells", fontsize_number=13, fontsize=13, number_format = "%.2f", breaks=seq(min, max, length.out=100),
               angle_col=0, legend=FALSE, silent=TRUE, color = col(100))
p2 <- pheatmap(t(res2), display_numbers=TRUE, number_color="black", border_color=NA, cluster_cols=FALSE, cluster_rows=FALSE, 
               main="Ductal Cells", fontsize_number=13, fontsize=13, number_format = "%.2f", breaks=seq(min, max, length.out=100),
               angle_col=0, legend=FALSE, silent=TRUE, color = col(100))

plot_list[[1]] <- p1[[4]]
plot_list[[2]] <- p2[[4]]


pdf("Final/CorrelationMatrix_OrganoidMarkerGenes_DownSample_Average100_Simple.pdf", width = 5, height = 3)
grid.arrange(arrangeGrob(grobs=plot_list[1:2], nrow=2))
dev.off()
