setwd('Z:/Van/20210531 - FD-seq technical replicates/data_analysis')

# Libraries
library(monocle)
library(ggplot2)

# Import data filtered by species from python
df1 <- read.table('Z:/Van/20210531 - FD-seq technical replicates/count_matrix/cDNA_rep1_dge.txt.gz', header=T, row.names=1)
df2 <- read.table('Z:/Van/20210531 - FD-seq technical replicates/count_matrix/cDNA_rep2_dge.txt.gz', header=T, row.names=1)

# Histogram
hist(colSums(df1))
hist(colSums(df1>0))
hist(colSums(df2))
hist(colSums(df2>0))

#*********************************************#
# # Discard cells with more than 5000 genes
# df1 <- df1[, colSums(df1>0)<=5e3]
# df2 <- df2[, colSums(df2>0)<=5e3]

# Discard cells with more than 15000 UMIs
df1 <- df1[, colSums(df1)<=15e3]
df2 <- df2[, colSums(df2)<=15e3]

# Discard cells with fewer than 1000 UMIs
df1 <- df1[, colSums(df1)>=1e3]
df2 <- df2[, colSums(df2)>=1e3]

# Keep genes detected in at least 5 cells
df1 <- df1[rowSums(df1>0)>=5,]
df2 <- df2[rowSums(df2>0)>=5,]

# Combine human data frames
df_all <- merge(df1, df2, by='row.names', all=TRUE)
row.names(df_all) <- df_all$Row.names
df_all <- df_all[,-1] # drop column Row.names
df_all[is.na(df_all)] <- 0 # fill missing values with 0

# Convert to CellDataSet object
# Phenotype data
pd1 <- data.frame(colnames(df1), row.names = colnames(df1))
pd2 <- data.frame(colnames(df2), row.names = colnames(df2))
colnames(pd1) <- 'cell.barcode'
colnames(pd2) <- 'cell.barcode'
# Annotate sample type
pd1$sample <- 'rep1'
pd2$sample <- 'rep2'
pd <- new('AnnotatedDataFrame', data=rbind(pd1,pd2))
# Feature data
fd <- data.frame(rownames(df_all), row.names = rownames(df_all))
colnames(fd) <- 'gene_short_name' # column name must be gene_short_name
fd <- new('AnnotatedDataFrame', data=fd)
cds_human <- newCellDataSet(Matrix(as.matrix(df_all),sparse=TRUE), phenoData=pd, featureData=fd, lowerDetectionLimit = 0.1, expressionFamily=negbinomial.size())

# Estimate size factors and dispersions
cds_human <- estimateSizeFactors(cds_human)
cds_human <- estimateDispersions(cds_human)

# # Export size factors
# write.table(sizeFactors(cds_human), file='Z:/Van/20190212 - Species mix live and fixed cells/all_human_sizefactors.txt', sep='\t', col.names=c('size_factor'))

# Extract normalized gene expression level
human_norm_gene_expr <- t(t(exprs(cds_human))/sizeFactors(cds_human))
write.table(as.matrix(human_norm_gene_expr), file=gzfile("Z:/Van/20210531 - FD-seq technical replicates/data_analysis/all_norm_gene_expr.txt.gz"), sep='\t')
