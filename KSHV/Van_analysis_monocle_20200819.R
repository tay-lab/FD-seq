setwd('Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/DGE files')

# Install bioconductor
# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite('monocle') # Install monocle

# Libraries
library(data.table) # for fast data import
library(Matrix) # sparse matrix
library(ggplot2) # for plotting
library(cowplot) # default theme for ggplot2
library(colorRamps) # color palettes
library(plyr)
library(dplyr) # for %% operation
# library(reticulate) # for monocle 3
library(monocle)
library(MASS) # negative binomial regression

# GO analysis (these 2 packages interfere with dplyr) -> use http://www.geneontology.org/
# library(GO.db)
# library(GOstats)

#library(limma)
#library(edgeR)

# Load custom fuctions
source('plot_pseudotime_heatmap_Van.R') # custom plot_pseudotime_heatmap, with labels_row = column name of fData to be used as row labels
source('plot_genes_violin_Van.R') # custom plot_genes_violin, with added option for fill color
source('plot_genes_in_pseudotime_Van.R') # custom plot_genes_in_pseudotime, with customized alpha and size

# Import data
raw_pos <- read.table('Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/DGE files/positive.dge.txt', row.names=1, header=TRUE)
raw_neg <- read.table('Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/DGE files/negative.dge.txt', row.names=1, header=TRUE)
# raw_pos <- fread('Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/DGE files/positive.dge.txt')

# Import list of KSHV genes
kshv_gtf <- fread('Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/DGE files/KSHV.gtf')
kshv_gtf <- kshv_gtf[kshv_gtf[,3][[1]]=='gene',] # keep the gene annotation only
kshv_gtf <- kshv_gtf[,9][[1]] # keep the last column only
kshv_gtf <- read.table(text=kshv_gtf, sep=';', stringsAsFactors=FALSE) # convert to columns
kshv_genes <- list(kshv_gtf[,6]) # Extract gene names
kshv_genes <- lapply(kshv_genes, FUN=function(x) gsub(' gene_name ','',x)) # remove common characters

# Add a prefix 'KSHV:' to viral genes
for (i in 1:length(kshv_genes[[1]])){
  if (any(rownames(raw_pos) == kshv_genes[[1]][i])){
    rownames(raw_pos)[rownames(raw_pos) == kshv_genes[[1]][i]] <- paste0('KSHV:', kshv_genes[[1]][i])
  }
  if (any(rownames(raw_neg) == kshv_genes[[1]][i])){
    rownames(raw_neg)[rownames(raw_neg) == kshv_genes[[1]][i]] <- paste0('KSHV:', kshv_genes[[1]][i])
  }
}

# KSHV genes by stage (ref:https://www.nature.com/articles/nri2076/figures/1 and https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1003847 )
# From 2nd ref: genes expressed 8 hours are IE, 24 or 24-48 hour are considered early, 48 hours are late
kshv_genes_latent <- list(paste0('KSHV:', c('K1','K2','K12','ORF71','ORF72','ORF73')))
kshv_genes_imme <- list(paste0('KSHV:', c('ORF11','ORF70','K4','K4.1','K4.2','K5','K6','ORF16','ORF45',
                                          'ORF50','K8','ORF57')))
kshv_genes_early <- list(paste0('KSHV:', c('ORF4','ORF6','K3','ORF17.5','ORF18','ORF34','ORF35','ORF36',
                                           'ORF37','ORF38','ORF39','ORF46','ORF47','ORF58','ORF59',
                                           'ORF60','ORF61','K14','ORF74')))
kshv_genes_late <- list(paste0('KSHV:', c('ORF8','ORF9','ORF10','K3A','ORF17','ORF21','ORF22','ORF23',
                                          'ORF24','ORF25','ORF26','ORF27','ORF28','ORF29','ORF30',
                                          'ORF31','ORF32','ORF33','ORF40','ORF41','ORF42','ORF43',
                                          'ORF44','ORF45.1','K8.1','ORF52','ORF53','ORF54','ORF55',
                                          'ORF56','K9','K10','K10.5','K11','ORF62','ORF65','ORF66',
                                          'ORF67','ORF67.5','ORF68','ORF69','ORF75')))
# K5/K6 antisense are late, not included

# QC
hist(colSums(raw_pos), xlab='Total UMI', main='Positive sample: Total UMI per cell', breaks=30)
hist(colSums(raw_pos>0), xlab='Total genes', main='Positive sample: Total genes per cell', breaks=30)
hist(colSums(raw_neg), xlab='Total UMI', main='Negative sample: Total UMI per cell', breaks=30)
hist(colSums(raw_neg>0), xlab='Total genes', main='Negative sample: Total genes per cell', breaks=30)

# Discard cells with >= 3000 genes
raw_pos <- raw_pos[, colSums(raw_pos>0) < 3e3]
raw_neg <- raw_neg[, colSums(raw_neg>0) < 3e3]

# Keep cells with UMI >= 1,000 and <= 10,000
raw_pos <- raw_pos[, colSums(raw_pos) >= 1e3 & colSums(raw_pos) <= 1e4]
raw_neg <- raw_neg[, colSums(raw_neg) >= 1e3 & colSums(raw_neg) <= 1e4]

# Keep genes detected in at least 5 cells
raw_pos <- raw_pos[rowSums(raw_pos>0) >=5,]
raw_neg <- raw_neg[rowSums(raw_neg>0) >=5,]

# Double check
hist(colSums(raw_pos), xlab='Total UMI', main='Filtered positive sample: Total UMI per cell', breaks=30)
hist(colSums(raw_pos>0), xlab='Total genes', main='Filtered positive sample: Total genes per cell', breaks=30)
hist(colSums(raw_neg), xlab='Total UMI', main='Filtered negative sample: Total UMI per cell', breaks=30)
hist(colSums(raw_neg>0), xlab='Total genes', main='Filtered negative sample: Total genes per cell', breaks=30)

# Histogram of % viral transcripts
pos.percent.kshv <- Matrix::colSums(raw_pos[grep('KSHV:',rownames(raw_pos)),])/Matrix::colSums(raw_pos) # percentage of KSHV transcripts
neg.percent.kshv <- Matrix::colSums(raw_neg[grep('KSHV:',rownames(raw_neg)),])/Matrix::colSums(raw_neg) # percentage of KSHV transcripts
ggplot() +
  geom_histogram(aes(x=pos.percent.kshv, y=(..count..)/sum(..count..),fill='Positive population'), colour='blue', alpha=0.5, bins=30) +
  geom_histogram(aes(x=neg.percent.kshv, y=(..count..)/sum(..count..),fill='Negative population'), colour='black', alpha=0.25, bins=30) +
  scale_fill_manual(name=NULL,values=c('Positive population'='blue','Negative population'='black')) +
  # scale_alpha_manual(values=c(0.5,0.25), guide=FALSE)
  # guides(colour=guide_legend((overisde.aes=list(alpha=c(0.5,0.25))))) +
  xlab('Proportion of KSHV transcripts') + ylab('Frequency') + ylim(0,0.8) +
  theme(legend.position=c(0.5,0.8)) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# Check if there's any common cell barcode
any(colnames(raw_pos) %in% colnames(raw_neg))

# Combine the two data sets
raw_combined <- merge(raw_pos, raw_neg, by='row.names', all=TRUE)
row.names(raw_combined) <- raw_combined$Row.names
raw_combined <- raw_combined[,-1] # drop column Row.names
raw_combined[is.na(raw_combined)] <- 0 # missing values should be 0 instead of NA

# Convert to CellDataSet object
# Phenotype data
pd <- data.frame(colnames(raw_combined), row.names = colnames(raw_combined))
colnames(pd) <- 'cell.barcode'
# pd$number.genes <- Matrix::colSums(raw_pos>0) # number of genes
pd$num_host_genes_expressed <- colSums(raw_combined[!grepl('KSHV:', row.names(raw_combined)),]>0) # number of detected host genes
pd$total_umi_host_genes <- colSums(raw_combined[!grepl('KSHV:', row.names(raw_combined)),]) # total UMI of host genes
pd$pct.kshv <- Matrix::colSums(raw_combined[grep('KSHV:',rownames(raw_combined)),])/Matrix::colSums(raw_combined) # % of KSHV transcripts
# Annotate sample type
pd$sample <- 'Positive'
pd$sample[(ncol(raw_pos)+1):ncol(raw_combined)] <- 'Negative'
pd <- new('AnnotatedDataFrame', data=pd)
# Feature data
fd <- data.frame(rownames(raw_combined), row.names = rownames(raw_combined))
colnames(fd) <- 'gene_short_name' # column name must be gene_short_name
fd$gene_type <- 'host' # host genes
fd$gene_type[grep('KSHV:', fd$gene_short_name)] <- 'KSHV:unannotated' # viral genes
# Annotate viral gene state
for (i in kshv_genes_imme[[1]]){
  if (i %in% fd$gene_short_name){
    fd$gene_type[fd$gene_short_name == i] <- 'KSHV:immediate'
  }
}
for (i in kshv_genes_early[[1]]){
  if (i %in% fd$gene_short_name){
    fd$gene_type[fd$gene_short_name == i] <- 'KSHV:early'
  }
}
for (i in kshv_genes_late[[1]]){
  if (i %in% fd$gene_short_name){
    fd$gene_type[fd$gene_short_name == i] <- 'KSHV:late'
  }
}
for (i in kshv_genes_latent[[1]]){
  if (i %in% fd$gene_short_name){
    fd$gene_type[fd$gene_short_name == i] <- 'KSHV:latent'
  }
}
# fd$num_cells_expressed <- Matrix::rowSums(raw_pos>0) # number of cells expressing the gene
fd <- new('AnnotatedDataFrame', data=fd)
cds.combined <- newCellDataSet(Matrix(as.matrix(raw_combined),sparse=TRUE), phenoData=pd, featureData=fd, lowerDetectionLimit = 0.1, expressionFamily=negbinomial.size())

# Estimate size factors and dispersions
cds.combined <- estimateSizeFactors(cds.combined)
cds.combined <- estimateDispersions(cds.combined)

# Statistics on number of genes/cell and number of cells/gene
cds.combined <- detectGenes(cds.combined, min_expr=0.1)

# use genes detected in at least 5% of the cells for clustering
fData(cds.combined)$use_for_ordering <- fData(cds.combined)$num_cells_expressed >= 0.05*ncol(cds.combined)
table(fData(cds.combined)$use_for_ordering) # summary

plot_pc_variance_explained(cds.combined, return_all=FALSE) # Find signifcant PCs
cds.combined <- reduceDimension(cds.combined, max_components=2, norm_method='log', num_dim=4, reduction_method='tSNE', verbose=TRUE) # Reduce dimension
cds.combined <- clusterCells(cds.combined, verbose=FALSE) # Cluster cells with default option
plot_cell_clusters(cds.combined)

# set thresholds for rho and delta to identify outliers at the top right corner, and # outliers = # clusters
plot_rho_delta(cds.combined, rho_threshold=17, delta_threshold=11) # 7 clusters, faster DDRTree dimension reduction for some reason
cds.combined <- clusterCells(cds.combined, rho_threshold=17, delta_threshold=11, skip_rho_sigma=TRUE)
plot_cell_clusters(cds.combined, color_by='Cluster')
plot_cell_clusters(cds.combined, color_by='sample')

# Export data for plotting in python
tsne_xy <- as.data.frame(t(reducedDimA(cds.combined)))
colnames(tsne_xy) <- c('tSNE1','tSNE2')
tsne_xy$sample <- pData(cds.combined)$sample
tsne_xy$pct_kshv <- pData(cds.combined)$pct.kshv # % viral transcripts
tsne_xy$K8.1 <- as.matrix(Matrix::t(exprs(cds.combined['KSHV:K8.1',]))/sizeFactors(cds.combined)) # normalized K8.1 expression
tsne_xy$pct_latent <- colSums(raw_combined[fData(cds.combined)$gene_type=='KSHV:latent',])/colSums(raw_combined) # % of latent viral transcripts
tsne_xy$pct_early <- colSums(raw_combined[fData(cds.combined)$gene_type=='KSHV:early',])/colSums(raw_combined) # % of early viral transcripts
tsne_xy$pct_imme <- colSums(raw_combined[fData(cds.combined)$gene_type=='KSHV:immediate',])/colSums(raw_combined) # % of immediate-early viral transcripts
tsne_xy$pct_late <- colSums(raw_combined[fData(cds.combined)$gene_type=='KSHV:late',])/colSums(raw_combined) # % of late viral transcripts

# How much of the viral transcripts in each cell is latent, immediate-early, early, and late
tsne_xy$prop_latent <- colSums(raw_combined[fData(cds.combined)$gene_type=='KSHV:latent',])/colSums(raw_combined[grepl('KSHV:',rownames(raw_combined)),])
tsne_xy$prop_early <- colSums(raw_combined[fData(cds.combined)$gene_type=='KSHV:early',])/colSums(raw_combined[grepl('KSHV:',rownames(raw_combined)),])
tsne_xy$prop_imme <- colSums(raw_combined[fData(cds.combined)$gene_type=='KSHV:immediate',])/colSums(raw_combined[grepl('KSHV:',rownames(raw_combined)),])
tsne_xy$prop_late <- colSums(raw_combined[fData(cds.combined)$gene_type=='KSHV:late',])/colSums(raw_combined[grepl('KSHV:',rownames(raw_combined)),])
tsne_xy$prop_unknown <- colSums(raw_combined[fData(cds.combined)$gene_type=='KSHV:unknown',])/colSums(raw_combined[grepl('KSHV:',rownames(raw_combined)),])

# Relative expression level of viral transcripts by kinetics
tsne_xy$rel_latent <- as.matrix(rowMeans(Matrix::t(exprs(cds.combined[fData(cds.combined)$gene_type=='KSHV:latent',]))/sizeFactors(cds.combined)))
tsne_xy$rel_early <- as.matrix(rowMeans(Matrix::t(exprs(cds.combined[fData(cds.combined)$gene_type=='KSHV:early',]))/sizeFactors(cds.combined)))
tsne_xy$rel_imme <- as.matrix(rowMeans(Matrix::t(exprs(cds.combined[fData(cds.combined)$gene_type=='KSHV:immediate',]))/sizeFactors(cds.combined)))
tsne_xy$rel_late <- as.matrix(rowMeans(Matrix::t(exprs(cds.combined[fData(cds.combined)$gene_type=='KSHV:late',]))/sizeFactors(cds.combined)))

write.csv(tsne_xy, file='tsne_by_sample.csv')

# # Plot change in early viral genes' expression level
# ggplot(tsne_xy, aes(x=pct_kshv*100,y=rel_imme)) +
#   geom_point(aes(color=sample), alpha=0.5, stroke=0) +
#   geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se=FALSE, colour='black', size=1.25) +
#   scale_color_brewer(palette='Dark2', direction=-1) +
#   xlab('Percentage of viral transcripts (%)') +
#   ylab('Normalized expression of early viral genes') +
#   ylim(-1,30) +
#   theme(legend.position = "none", axis.title=element_text(size=12))
# 
# # Plot change in % late
# ggplot(tsne_xy, aes(x=pct_kshv*100,y=pct_late*100)) +
#   geom_point(aes(color=sample), alpha=0.5, stroke=0, size=1.3) +
#   # geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se=FALSE, colour='black', size=1.25) +
#   scale_color_brewer(palette='Dark2', direction=-1) +
#   xlab('Percentage of viral transcripts (%)') +
#   ylab('Percentage of late\nviral transcripts (%)') +
#   ylim(-1,30) +
#   theme_classic() +
#   theme(legend.position = "none", axis.title=element_text(size=12))
# ggsave('Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/figures/late_vs_pct_kshv.pdf', width=3,height=3,units='in')
# 
# # Plot change in % early
# ggplot(tsne_xy, aes(x=pct_kshv*100,y=pct_early*100)) +
#   geom_point(aes(color=sample), alpha=0.5, stroke=0, size=1.3) +
#   geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se=FALSE, colour='black', size=1.25) +
#   scale_color_brewer(palette='Dark2', direction=-1) +
#   xlab('Percentage of viral transcripts (%)') +
#   ylab('Percentage of early viral transcripts (%)') +
#   theme(legend.position = "none", axis.title=element_text(size=12))
# ggsave('Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/figures/early_vs_pct_kshv.pdf', width=3,height=3,units='in')
# 
# # Plot change in % immediate
# ggplot(tsne_xy, aes(x=pct_kshv*100,y=pct_imme*100)) +
#   geom_point(aes(color=sample), alpha=0.5, stroke=0, size=1.3) +
#   geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se=FALSE, colour='black', size=1.25) +
#   scale_color_brewer(palette='Dark2', direction=-1) +
#   xlab('Percentage of viral transcripts (%)') +
#   ylab('Percentage of immediate-early viral transcripts (%)') +
#   theme(legend.position = "none", axis.title=element_text(size=12))
# ggsave('Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/figures/immediate_vs_pct_kshv.pdf', width=3,height=3,units='in')
# 
# # Plot change in % latent
# ggplot(tsne_xy, aes(x=pct_kshv*100,y=pct_latent*100)) +
#   geom_point(aes(color=sample), alpha=0.5, stroke=0, size=1.3) +
#   geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se=FALSE, colour='black', size=1.25) +
#   scale_color_brewer(palette='Dark2', direction=-1) +
#   xlab('Percentage of viral transcripts (%)') +
#   ylab('Percentage of latent viral transcripts (%)') +
#   theme(legend.position = "none", axis.title=element_text(size=12))
# ggsave('Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/figures/latent_vs_pct_kshv.pdf', width=3,height=3,units='in')

# Plot of K8.1 expression level
plot_genes_violin(cds.combined['KSHV:K8.1',], grouping='sample', color_by='sample') +
  scale_fill_brewer(palette='Dark2',direction=-1) +
  xlab('') + theme(legend.position = "none", axis.title=element_text(size=12))
plot_genes_jitter(cds.combined['KSHV:K8.1',], grouping='sample', color_by='sample', cell_size=1.15) +
  scale_color_brewer(palette='Dark2',direction=-1) +
  xlab('') + theme(legend.position = "none", axis.title=element_text(size=12), axis.text=element_text(size=12,color="black"))
ggsave('Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/figures/k8.1_jitter.pdf', width=3,height=3,units='in')

# Differential expressed genes between samples
sample.DEgenes <- differentialGeneTest(cds.combined[row.names(cds.combined)[fData(cds.combined)$use_for_ordering],], fullModelFormulaStr='~sample', cores=1)
# sample.DEgenes <- sample.DEgenes[sample.DEgenes$use_for_ordering,] # only look at genes expressed in >= 5% of cells
sample.DEgenes %>% arrange(qval) %>% head() # summary

# Plot expression level between two populations
sample.DEgenes %>% arrange(qval) %>% head() %>% dplyr::select(gene_short_name) -> top.sample.DEgenes
top.sample.DEgenes <- as.character(top.sample.DEgenes$gene_short_name)
plot_genes_violin_Van(cds.combined[top.sample.DEgenes,], grouping='sample', ncol=length(top.sample.DEgenes)/2, nrow=2, fill_color='cornflowerblue')

# Top 6 host genes
subset(sample.DEgenes, !grepl('KSHV:', gene_short_name)) %>% arrange(qval) %>% head() %>% dplyr::select(gene_short_name) -> top.sample.DEhostgenes
top.sample.DEhostgenes <- as.character(top.sample.DEhostgenes$gene_short_name)
plot_genes_violin_Van(cds.combined[top.sample.DEhostgenes,], grouping='sample', ncol=length(top.sample.DEhostgenes)/2, nrow=2, fill_color='cornflowerblue')

# Top 6 viral genes
subset(sample.DEgenes, grepl('KSHV:', gene_short_name)) %>% arrange(qval) %>% head() %>% dplyr::select(gene_short_name) -> top.sample.DEviralgenes
top.sample.DEviralgenes <- as.character(top.sample.DEviralgenes$gene_short_name)
plot_genes_violin_Van(cds.combined[top.sample.DEviralgenes,], grouping='sample', ncol=length(top.sample.DEviralgenes)/2, nrow=2, fill_color='cornflowerblue')

# DE viral genes
sample.DEgenes <- differentialGeneTest(cds.combined[grepl('KSHV:',rownames(cds.combined)),], fullModelFormulaStr='~Cluster')
sample.DEgenes %>% arrange(qval) %>% head() # summary

# # Gene dispersion
# pos.gene.disp <- dispersionTable(cds.combined)
# 
# # Only use genes with mean expression >=0.1 for clustering
# pos.clustering.genes <- pos.gene.disp[pos.gene.disp$mean_expression >=0.1,]
# cds.combined <- setOrderingFilter(cds.combined, pos.clustering.genes$gene_id)
# plot_ordering_genes(cds.combined)
# 
# # Plot major PCs
# plot_pc_variance_explained(cds.combined, return_all=FALSE)
# 
# # Reduce dimensions using t-SNE and cluster
# cds.combined <- reduceDimension(cds.combined, max_components=2, num_dim=10, reduction_method='tSNE', verbose=TRUE)
# cds.combined <- clusterCells(cds.combined, num_clusters=5)
# plot_cell_clusters(cds.combined, x=1, y=2)

################
# Differential expression ~percent.kshv (regressing out num_host_genes_expressed doesn't seem to have an effect)
################
cds.combined.percent.kshv <- cds.combined

# Assign percent of viral genes as pseudotime
pData(cds.combined.percent.kshv)$Pseudotime <- pData(cds.combined.percent.kshv)$pct.kshv*100

# DE on genes expressed in at least 5% of the cells only
cds.combined.percent.kshv <- detectGenes(cds.combined.percent.kshv, min_expr=0.1)
genes.to.test <- subset(fData(cds.combined.percent.kshv), fData(cds.combined.percent.kshv)$num_cells_expressed >= 0.05*ncol(cds.combined.percent.kshv))
genes.to.test <- row.names(genes.to.test)
percent.kshv.DEgenes <- differentialGeneTest(cds.combined.percent.kshv[genes.to.test,], fullModelFormulaStr='~sm.ns(Pseudotime)', cores=detectCores()-2)

# Export DE genes to csv file
write.csv(percent.kshv.DEgenes, "DEgenes_percent_kshv.csv")

# percent.kshv.DEgenes <- differentialGeneTest(cds.combined.percent.kshv[genes.to.test,],
#                                              fullModelFormulaStr='~sm.ns(Pseudotime)+total_umi_host_genes', reducedModelFormulaStr='~total_umi_host_genes',
#                                              cores=detectCores()-2)
# regressing out anything related to host genes/number of genes is not as good as actual downsampling

# # Use top 1000 DE genes for ordering
# my.ordering.genes <- row.names(percent.kshv.DEgenes)[order(percent.kshv.DEgenes$qval)[1:1000]]
# cds.combined.percent.kshv <- setOrderingFilter(cds.combined.percent.kshv, ordering_genes=my.ordering.genes)
# cds.combined.percent.kshv <- reduceDimension(cds.combined.percent.kshv, method='DDRTree')
# cds.combined.percent.kshv <- orderCells(cds.combined.percent.kshv)
# plot_cell_trajectory(cds.combined.percent.kshv)

# # Keep significant genes
# percent.kshv.DEgenes <- subset(percent.kshv.DEgenes, qval < 0.05)

# Significant, host, non-mitochondria DE genes
host.DEgenes <- subset(percent.kshv.DEgenes, qval <= 0.05 & gene_type=='host' & !grepl('MT-', gene_short_name))
# host.DEgenes <- subset(percent.kshv.DEgenes, qval < 1e-2 & gene_type=='host')
host.DEgenes <- as.character(host.DEgenes$gene_short_name)
plot.new()
my.heatmap <- plot_pseudotime_heatmap(cds.combined.percent.kshv[host.DEgenes,], num_clusters=6, use_gene_short_name=TRUE, show_rownames=FALSE, return_heatmap=TRUE)

my.cluster <- cutree(my.heatmap$tree_row, 6)
my.up.genes <- my.cluster[my.cluster==6] # extract genes from clusters 4 and 6
# plot.new()
# plot_pseudotime_heatmap(cds.pos.percent.kshv[my.up.genes,], num_clusters=3, use_gene_short_name=TRUE, show_rownames=TRUE)

my.up.genes.arranged <- percent.kshv.DEgenes[names(my.up.genes),c('gene_short_name','qval',"num_cells_expressed")] %>% arrange(qval)
print(my.up.genes.arranged)
my.up.genes.arranged <- as.character(my.up.genes.arranged$gene_short_name)
pdf('genes_upregulated_in_combined_dataset.pdf', width=12)
for (i in my.up.genes.arranged){
  print(plot_genes_in_pseudotime(cds.combined.percent.kshv[i,], color_by='sample') + xlab('Percentage of viral transcripts'))
}
dev.off()

for (i in c('ISCU','CDH1','CORO1C','TMEM119')){
  print(plot_genes_in_pseudotime_Van(cds.combined.percent.kshv[i,], color_by='sample', cell_size=1.3) +
          xlab('Percentage of viral transcripts (%)') + ylab('Relative expression') +
          scale_color_brewer(palette='Dark2', direction=-1) +
          theme_classic() +
          theme(legend.position = "none", axis.title=element_text(size=12), axis.text=element_text(size=12,color="black")))
  ggsave(paste0('Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/figures/',i,'.pdf'), width=3,height=3,units='in')
}

# Extract scaled expression level of positive population
scaled_expr <- Matrix::t(Matrix::t(exprs(cds.combined.percent.kshv))/sizeFactors(cds.combined.percent.kshv))
scaled_expr_pos <- scaled_expr[,pData(cds.combined.percent.kshv)$sample=='Positive']

# Export scaled expression level of KSHV genes
scaled_expr_pos_kshv <- scaled_expr_pos[grepl('KSHV:',rownames(scaled_expr_pos)),]
write.csv(as.matrix(scaled_expr_pos_kshv), file='scaled_expr_pos_kshv.csv')
# Export % viral transcripts
write.csv(pData(cds.combined.percent.kshv)[pData(cds.combined.percent.kshv)$sample=='Positive','pct.kshv',drop=FALSE],
          file='pct_viral_pos_kshv.csv')

# Export scaled expression level of host genes that are expressed in >=5% of the positive cells
write.csv(as.matrix(scaled_expr[!grepl('KSHV:',rownames(scaled_expr_pos)) &
                                  (fData(cds.combined.percent.kshv)$num_cells_expressed >= 0.05*ncol(scaled_expr_pos))]),
          file='scaled_expr_pos_host.csv')


################
# Cluster positive populations that have high % of viral transcripts
################
raw_pos2 <- raw_pos[,pos.percent.kshv>=0.3]

# Convert to cell data set
# Phenotype data
pd <- data.frame(colnames(raw_pos2), row.names = colnames(raw_pos2))
colnames(pd) <- 'cell.barcode'
# pd$number.genes <- Matrix::colSums(raw_pos>0) # number of genes
pd$num_host_genes_expressed <- colSums(raw_pos2[!grepl('KSHV:', row.names(raw_pos2)),]>0) # number of detected host genes
pd$total_umi_host_genes <- colSums(raw_pos2[!grepl('KSHV:', row.names(raw_pos2)),]) # total UMI of host genes
pd$pct.kshv <- Matrix::colSums(raw_pos2[grep('KSHV:',rownames(raw_pos2)),])/Matrix::colSums(raw_pos2) # % of KSHV transcripts
# Annotate sample type
pd$sample <- 'Positive'
pd <- new('AnnotatedDataFrame', data=pd)
# Feature data
fd <- data.frame(rownames(raw_pos2), row.names = rownames(raw_pos2))
colnames(fd) <- 'gene_short_name' # column name must be gene_short_name
fd$gene_type <- 'host' # host genes
fd$gene_type[grep('KSHV:', fd$gene_short_name)] <- 'KSHV:unknown' # viral genes
# Annotate viral gene state
for (i in kshv_genes_early[[1]]){
  if (i %in% fd$gene_short_name){
    fd$gene_type[fd$gene_short_name == i] <- 'KSHV:early'
  }
}
for (i in kshv_genes_imme[[1]]){
  if (i %in% fd$gene_short_name){
    fd$gene_type[fd$gene_short_name == i] <- 'KSHV:immediate'
  }
}
for (i in kshv_genes_late[[1]]){
  if (i %in% fd$gene_short_name){
    fd$gene_type[fd$gene_short_name == i] <- 'KSHV:late'
  }
}
for (i in kshv_genes_latent[[1]]){
  if (i %in% fd$gene_short_name){
    fd$gene_type[fd$gene_short_name == i] <- 'KSHV:latent'
  }
}
# fd$num_cells_expressed <- Matrix::rowSums(raw_pos>0) # number of cells expressing the gene
fd <- new('AnnotatedDataFrame', data=fd)
cds.pos2 <- newCellDataSet(Matrix(as.matrix(raw_pos2),sparse=TRUE), phenoData=pd, featureData=fd, lowerDetectionLimit = 0.1, expressionFamily=negbinomial.size())

# Estimate size factors and dispersions
cds.pos2 <- estimateSizeFactors(cds.pos2)
cds.pos2 <- estimateDispersions(cds.pos2)

# Statistics on number of genes/cell and number of cells/gene
cds.pos2 <- detectGenes(cds.pos2, min_expr=0.1)

# use genes detected in at least 5% of the cells for clustering
fData(cds.pos2)$use_for_ordering <- fData(cds.pos2)$num_cells_expressed >= 0.05*ncol(cds.pos2)
table(fData(cds.pos2)$use_for_ordering) # summary

# Cluster cells with default option
plot_pc_variance_explained(cds.pos2, return_all=FALSE) # Find signifcant PCs
cds.pos2 <- reduceDimension(cds.pos2, max_components=2, norm_method='log', num_dim=5, reduction_method='tSNE', verbose=TRUE) # Reduce dimension
cds.pos2 <- clusterCells(cds.pos2, verbose=FALSE)
plot_cell_clusters(cds.pos2)

# set thresholds for rho and delta to identify outliers at the top right corner, and # outliers = # clusters
plot_rho_delta(cds.pos2, rho_threshold=14, delta_threshold=11) # 7 clusters, faster DDRTree dimension reduction for some reason
cds.pos2 <- clusterCells(cds.pos2, rho_threshold=14, delta_threshold=11, skip_rho_sigma=TRUE)
plot_cell_clusters(cds.pos2, color_by='Cluster')

# DE genes
pos2.DEgenes <- differentialGeneTest(cds.pos2[row.names(cds.pos2)[fData(cds.pos2)$use_for_ordering],], fullModelFormulaStr='~Cluster', cores=1)
pos2.DEgenes %>% arrange(qval) %>% head() # summary

# ######### Negative binomial regression on host genes of positive population
# # Extract host genes
# scaled_expr_pos_hostgenes <- scaled_expr_pos[!grepl('KSHV:',rownames(scaled_expr_pos)),]
# # Keep cells detected in at least 5% of the cells
# scaled_expr_pos_hostgenes <- scaled_expr_pos_hostgenes[Matrix::rowSums(scaled_expr_pos_hostgenes>0)>=0.05*ncol(scaled_expr_pos_hostgenes),]
# # Round the count data
# scaled_expr_pos_hostgenes <- round(scaled_expr_pos_hostgenes)
# 
# # Keep cells with >= 20% of viral transcripts
# pos_pct_kshv <- pData(cds.combined.percent.kshv)[pData(cds.combined.percent.kshv)$sample=='Positive','pct.kshv']
# scaled_expr_pos_hostgenes <- scaled_expr_pos_hostgenes[,pos_pct_kshv>=0.2]
# 
# # Negative binomial regression
# my_nbr <- data.frame(p_value=numeric(nrow(scaled_expr_pos_hostgenes))+1,
#                      slope=numeric(nrow(scaled_expr_pos_hostgenes)))
# rownames(my_nbr) <- rownames(scaled_expr_pos_hostgenes)
# nbr_pct_kshv <- pos_pct_kshv[pos_pct_kshv>=0.2]
# for (i in rownames(my_nbr)){
#   my_glm <- glm.nb(scaled_expr_pos_hostgenes[i,]~nbr_pct_kshv)
#   my_nbr[i,'p_value'] <- summary(my_glm)$coef[2,4]
#   my_nbr[i,'slope'] <- summary(my_glm)$coef[2,1]
# }
