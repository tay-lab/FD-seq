#---- Libraries
library(data.table)
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Seurat)
library(viridis)

# ggplot theme
my.theme <- theme(axis.title = element_text(size = 12))

# Change directory
setwd('Z:/Van/20200526 - Drop-seq fixed A549 OC43/data_analysis')

# Data analysis pipeline follows exactly: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

# Import mock data
df <- read.table("Z:/Van/20200526 - Drop-seq fixed A549 OC43/alignment_result/mock_dge.txt.gz", sep='\t', row.names = 1, header = TRUE)
# Filter data
# Keep cells with more than 1,000 UMIs, less than 20,000 UMIs
df <- df[,(colSums(df)>=1000) & (colSums(df)<=20000)]
# Keep genes detected in at least 5 cells
df <- df[rowSums(df>0) >= 5,]
# Convert into sparse seurat object
mock <- CreateSeuratObject(as(as.matrix(df), "sparseMatrix"))
mock$sample <- 'mock'

# MOI01
df <- read.table("Z:/Van/20200526 - Drop-seq fixed A549 OC43/alignment_result/MOI01_dge.txt.gz", sep='\t', row.names = 1, header = TRUE)
df <- df[,(colSums(df)>=1000) & (colSums(df)<=20000)]
df <- df[rowSums(df>0) >= 5,]
MOI01 <- CreateSeuratObject(as(as.matrix(df), "sparseMatrix"))
MOI01$sample <- 'MOI 0.1'

# MOI1
df <- read.table("Z:/Van/20200526 - Drop-seq fixed A549 OC43/alignment_result/MOI1_dge.txt.gz", sep='\t', row.names = 1, header = TRUE)
df <- df[,(colSums(df)>=1000) & (colSums(df)<=20000)]
df <- df[rowSums(df>0) >= 5,]
MOI1 <- CreateSeuratObject(as(as.matrix(df), "sparseMatrix"))
MOI1$sample <- 'MOI 1'

# a549 the 3 data sets
# a549 <- merge(mock, y=c(MOI01, MOI1), add.cell.ids = c('mock','MOI0.1','MOI1'))
a549 <- merge(mock, y=MOI1, add.cell.ids = c('mock','MOI1'))

# Remove objects to free up memory
rm(df, mock, MOI01, MOI1)

# Percent mitochondria
a549$pct.mt <- PercentageFeatureSet(a549, pattern = "^MT-")
# Filter cells with more than 20% mitochondria
a549 <- subset(a549, pct.mt < 20)

# Plot number of genes and UMIs
VlnPlot(a549, features = c('nCount_RNA','nFeature_RNA'))

# Count data
a549.counts <- GetAssayData(a549, slot = 'counts')

# Get the list of OC43 viral genes
viral.genes <- grep('OC43:', rownames(a549), value = TRUE)

# Export viral genes, and use python to plot % of OC43-positive cells
write.csv(a549.counts[viral.genes,], file = "all_samples_viral_genes.csv")

# Calculate percentage of viral genes
a549$pct.viral <- colSums(a549.counts[viral.genes,])/colSums(a549.counts)*100
# Log of percentage of viral genes
a549$pct.viral.log <- log10(a549$pct.viral + 1)
# Calculate total UMI of viral genes
a549$total.viral <- colSums(a549.counts[viral.genes,])

# Normalize, scale, cluster data
a549 <- NormalizeData(a549)
a549 <- FindVariableFeatures(a549, selection.method = 'vst')
a549 <- ScaleData(a549)
a549 <- RunPCA(a549, features = VariableFeatures(a549))
# DimPlot(a549, reduction = "pca", label = TRUE, ident = 'sample')
# DimHeatmap(a549, dims = 1, cells = 500, balanced = TRUE)
# ElbowPlot(a549)
a549 <- FindNeighbors(a549, dims = 1:10)
a549 <- FindClusters(a549, resolution = 0.2)
a549 <- RunTSNE(a549, dims = 1:10)
DimPlot(a549, reduction = 'tsne') + xlab("t-SNE 1") + ylab("t-SNE 2") + my.theme
# DimPlot(a549, reduction = 'tsne', group.by = 'sample')
# # a549 <- RunUMAP(a549, dims = 1:10)
# # DimPlot(a549, reduction = 'umap')
# # DimPlot(a549, reduction = 'umap', group.by = 'sample')
# FeaturePlot(a549, features = c('pct.viral','total.viral'), reduction = 'tsne')
# VlnPlot(a549, features = c('pct.viral','total.viral'))
# VlnPlot(a549, features = c('pct.viral','total.viral'), group.by = 'sample')
# 
# # Cluster markers
# cluster.markers <- FindAllMarkers(a549, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# # Retain significant genes only
# cluster.markers <- cluster.markers[cluster.markers$p_val_adj<0.05,]
# # write.csv(cluster.markers, file = "5_clusters_markers.csv")
# cluster.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# Cell cycle scoring
a549 <- CellCycleScoring(a549, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
a549$cc.ident <- Idents(a549)
DimPlot(a549, reduction = 'tsne') + xlab("t-SNE 1") + ylab("t-SNE 2") + my.theme
ggsave("Z:/Van/20200526 - Drop-seq fixed A549 OC43/data_analysis/tsne_cellcycle.png", dpi = 600, width = 3.8, height = 3, units = 'in')


# # Regress out cell cycle
# a549$cc.difference <- a549$S.Score - a549$G2M.Score
# a549 <- ScaleData(a549, vars.to.regress = "cc.difference")
# a549 <- NormalizeData(a549)
# a549 <- FindVariableFeatures(a549, selection.method = 'vst')
a549 <- ScaleData(a549, vars.to.regress = c("S.Score","G2M.Score"))
a549 <- RunPCA(a549, features = VariableFeatures(a549))
ElbowPlot(a549)
a549 <- FindNeighbors(a549, dims = 1:10)
a549 <- FindClusters(a549, resolution = 0.2)
a549 <- RunTSNE(a549, dims = 1:10)
# VlnPlot(a549, features = c('nCount_RNA','nFeature_RNA'))

# Rename cluster
Idents(a549) <- plyr::mapvalues(Idents(a549), from = c(0,1,2), to = c(1,0,2))
Idents(a549) <- factor(Idents(a549), levels = c("0","1","2"))

# t-SNE plots
DimPlot(a549, reduction = 'tsne', group.by = 'cc.ident') + xlab("t-SNE 1") + ylab("t-SNE 2") + my.theme
ggsave("Z:/Van/20200526 - Drop-seq fixed A549 OC43/data_analysis/tsne_no_cellcycle.png", dpi = 600, width = 3.8, height = 3, units = 'in')
DimPlot(a549, reduction = 'tsne') + xlab("t-SNE 1") + ylab("t-SNE 2") + my.theme
ggsave("Z:/Van/20200526 - Drop-seq fixed A549 OC43/data_analysis/tsne_clusters.png", dpi = 600, width = 3.8, height = 3, units = 'in')
DimPlot(a549, reduction = 'tsne', group.by = 'sample') + xlab("t-SNE 1") + ylab("t-SNE 2") + my.theme
ggsave("Z:/Van/20200526 - Drop-seq fixed A549 OC43/data_analysis/tsne_sample.png", dpi= 600, width = 4.2, height = 3, units = 'in')

# Plot percentage of viral gene
VlnPlot(a549, features = "pct.viral") + xlab("Cluster") + ylab("Percentage of viral transcripts") + my.theme
ggsave("Z:/Van/20200526 - Drop-seq fixed A549 OC43/data_analysis/violin_pct_viral.png", dpi = 600, width = 3.5, height = 3.2, units = 'in')
VlnPlot(a549, features = "pct.viral.log") + xlab("Cluster") + ylab("log10 (% viral transcripts + 1)") + my.theme
ggsave("Z:/Van/20200526 - Drop-seq fixed A549 OC43/data_analysis/violin_pct_viral_log.png", dpi = 600, width = 3.5, height = 3.2, units = 'in')
# FeaturePlot(a549, features = c('pct.viral'), reduction = 'tsne') + xlab("t-SNE 1") + ylab("t-SNE 2") + my.theme
# ggsave("Z:/Van/20200526 - Drop-seq fixed A549 OC43/data_analysis/tsne_pct_viral.png", dpi = 600, width = 3.8, height = 3, units = 'in')

# Plot total viral gene count
a549$total.viral.log <- log10(a549$total.viral+1)
VlnPlot(a549, features=c("total.viral.log"), group.by = 'sample') + xlab("Sample") + ylab("log10 (total viral UMI + 1)") + my.theme
ggsave("Z:/Van/20200526 - Drop-seq fixed A549 OC43/data_analysis/violin_total_viral_log.png", dpi = 600, width = 3, height = 3.5, units = 'in')
FeaturePlot(a549, features = c("pct.viral.log","total.viral.log"))
ggsave("Z:/Van/20200526 - Drop-seq fixed A549 OC43/data_analysis/tsne_pct_viral_log.png", dpi = 600, width = 7.5, height = 3.5, units = 'in')

# Get cluster marker
cluster.markers.cc <- FindAllMarkers(a549, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers.cc <- cluster.markers.cc[cluster.markers.cc$p_val_adj<0.05,]
# write.csv(cluster.markers.cc, file = "3_clusters_markers_cellcycle.csv")
cluster.markers.cc %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# Infected cells
infected <- subset(a549, cells = colnames(a549)[Idents(a549)==2])
VlnPlot(infected, features = c('pct.viral', 'total.viral'))
# VlnPlot(infected, features = viral.genes, nrow = 2)

# Export the relative expression of viral genes in cluster 2 to plot in python
infected.data <- GetAssayData(infected, slot = 'data')
infected.data <- infected.data[viral.genes,]
infected.data <- as.data.frame(infected.data)
infected.data['pct.viral',] <- infected$pct.viral[colnames(infected.data)]
write.csv(t(infected.data), file = "infected_relative_viral_genes.csv")

# Export the relative expression of viral genes in MOI 1 sample to plot in python
moi1.data <- GetAssayData(subset(a549, subset = (sample == 'MOI 1')), slot = 'data')
moi1.data <- as.data.frame(t(moi1.data[viral.genes,]))
moi1.data$cluster <- Idents(a549)[rownames(moi1.data)] # add info about cluster
moi1.data$total.viral.log <- a549$total.viral.log[rownames(moi1.data)]
write.csv(moi1.data, file = "moi1_relative_viral_genes.csv")
# High subpopulation
10^median(moi1.data$total.viral.log[moi1.data$total.viral.log>2])-1
# Low subpopulation
10^median(moi1.data$total.viral.log[(moi1.data$total.viral.log<2)&(moi1.data$total.viral.log>0)])-1

# # Find genes that are correlated with total.viral.log
# moi1.data <- GetAssayData(subset(a549, subset = (sample == 'MOI 1')), slot = 'data')
# my.lm.results <- data.frame(slope=numeric(nrow(moi1.data)),
#                     pval=numeric(nrow(moi1.data)) + 1,
#                     row.names = rownames(moi1.data))
# temp.total.viral.log <- a549$total.viral.log[colnames(moi1.data)]
# for (gene_i in rownames(moi1.data))
# {
#   temp.lm <- lm(as.numeric(moi1.data[gene_i,])~as.numeric(temp.total.viral.log))
#   my.lm.results[gene_i,'slope'] <- coef(temp.lm)[2]
#   my.lm.results[gene_i,'pval'] <- summary(temp.lm)$coefficients[2,4]
# }
# # Replace pval==NA by 1
# my.lm.results[is.na(my.lm.results$pval),'pval'] <- 1
# my.lm.results$pval.adj <- p.adjust(my.lm.results$pval, method = 'fdr')
# my.lm.results[my.lm.results$pval.adj<1e-150,'pval.adj'] <- 1e-150
# plot(my.lm.results$slope, -log10(my.lm.results$pval.adj))
# pos.only <- my.lm.results[(my.lm.results$slope > 0.25) & (my.lm.results$pval.adj<0.05),]
# # write.csv(rownames(pos.only), file = "positive_lm.csv")

# Plot 4 markers of cluster 0
# FeaturePlot(a549, features = c("CXCL1","CXCL5","NFKBIA","CCL2"))
VlnPlot(a549, features = c("CXCL1","CXCL5"), ncol = 2, pt.size = 0.1)
ggsave("Z:/Van/20200526 - Drop-seq fixed A549 OC43/data_analysis/immune_genes1_violin.png", dpi = 600, width = 6.2, height = 3.2, units = 'in')
VlnPlot(a549, features = c("CCL2","NFKBIA"), ncol = 2, pt.size = 0.1)
ggsave("Z:/Van/20200526 - Drop-seq fixed A549 OC43/data_analysis/immune_genes2_violin.png", dpi = 600, width = 6.2, height = 3.2, units = 'in')

# Plot all signaling pathway genes on t-sne
FeaturePlot(a549, features = c("CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL8", "NFKBIA", "CCL2","CEBPB","PTGS2"), ncol=3)
ggsave("Z:/Van/20200526 - Drop-seq fixed A549 OC43/data_analysis/tsne_cluster1_9markers.png", dpi = 600, width = 9, height = 8, units = 'in')
# VlnPlot(a549, features = c("CXCL1","CXCL5","NFKBIA","CCL2"), ncol = 2, group.by = 'sample', pt.size = 0)

#----- Pseudotime
# library(monocle)
# 
# # Obtain MOI 1 cells
# moi1.counts <- GetAssayData(subset(a549, subset = (sample == 'MOI 1')), slot = 'counts')
# 
# # Phenotype data
# pd <- data.frame(colnames(moi1.counts))
# rownames(pd) <- colnames(moi1.counts)
# colnames(pd) <- 'cell.barcode'
# pd <- new('AnnotatedDataFrame', data=pd)
# # Feature data
# fd <- data.frame(rownames(moi1.counts))
# rownames(fd) <- rownames(moi1.counts)
# colnames(fd) <- 'gene_short_name'
# fd$gene_type <- 'host'
# fd[grep('OC43:',rownames(fd)),'gene_type'] <- 'OC43'
# fd <- new('AnnotatedDataFrame', data=fd)
# # CellDataSet
# cds.moi1 <- newCellDataSet(moi1.counts, phenoData=pd, featureData=fd, lowerDetectionLimit = 0.1, expressionFamily=negbinomial.size())
# 
# # Estimate size factors and dispersions
# cds.moi1 <- estimateSizeFactors(cds.moi1)
# cds.moi1 <- estimateDispersions(cds.moi1)
# 
# # Perform pseudotime analysis
# cds.moi1 <- detectGenes(cds.moi1, min_expr=0.1)
# fData(cds.moi1)$use_for_ordering <- fData(cds.moi1)$num_cells_expressed >= 0.05*ncol(cds.moi1)
# plot_pc_variance_explained(cds.moi1, return_all=FALSE) # Find signifcant PCs
# cds.moi1 <- reduceDimension(cds.moi1, max_components=2, norm_method='log', num_dim=5, reduction_method='tSNE', verbose=TRUE) # Reduce dimension
# cds.moi1 <- clusterCells(cds.moi1, verbose=FALSE) # Cluster cells with default option
# plot_cell_clusters(cds.moi1)
# plot_rho_delta(cds.moi1, rho_threshold=17, delta_threshold=15) # 4 clusters
# cds.moi1 <- clusterCells(cds.moi1, rho_threshold=17, delta_threshold=15, skip_rho_sigma=TRUE)
# plot_cell_clusters(cds.moi1)
# # DE genes
# cds.moi1.de <- differentialGeneTest(cds.moi1, fullModelFormulaStr = "~Cluster", cores = 2)
# cds.top1000.de <- rownames(cds.moi1.de)[order(cds.moi1.de$qval)][1:1000]
# cds.moi1 <- setOrderingFilter(cds.moi1, ordering_genes = cds.top1000.de)
# cds.moi1 <- reduceDimension(cds.moi1, method = "DDRTree")
# cds.moi1 <- orderCells(cds.moi1)
# plot_cell_trajectory(cds.moi1, color_by = "Cluster")
# 
# # Genes that change with pseudotime
# cds.moi1.pseudotime.de <- differentialGeneTest(cds.moi1, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 2)
# # plot_genes_in_pseudotime(cds.moi1[viral.genes,], color_by = "Cluster")
# 
# # Add abundance of viral transcripts
# pData(cds.moi1)$total.viral.log <- log10(colSums(moi1.counts[grep('OC43:',rownames(moi1.counts)),])+1)
# plot_cell_trajectory(cds.moi1, color_by = "total.viral.log")
# 
# # Analyze branch 1
# BEAM_1 <- BEAM(cds.moi1, branch_point = 1, cores = 2)
# BEAM_1 <- BEAM_1[order(BEAM_1$qval),]
# plot_genes_branched_heatmap(cds.moi1[rownames(head(BEAM_1,20)),],
#                             branch_point = 1, num_clusters = 4, use_gene_short_name = TRUE, show_rownames = TRUE)
# BEAM_2 <- BEAM(cds.moi1, branch_point = 2, cores = 2)
# BEAM_2 <- BEAM_2[order(BEAM_2$qval),]
# plot_genes_branched_heatmap(cds.moi1[rownames(head(BEAM_2,20)),],
#                             branch_point = 1, num_clusters = 4, use_gene_short_name = TRUE, show_rownames = TRUE)