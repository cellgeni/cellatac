#### Script to identify clusters of cells for a per cluster-based peak calling
#
# Adapted from pinellolab/scATAC-benchmarking Github repo
#   > Extra/10x_PBMC_5k/test_pseudobulk/Cusanovich2018_10xpbmc5k_idenfify_clades.ipynb 
# Based on http://atlas.gs.washington.edu/fly-atac/docs/
#


# TODO: make this useful also outside NF.
# This means supplying it perhaps with command line arguments.


theargs <- R.utils::commandArgs(asValues=TRUE)

N_clades = 10 # Number of expected cell clusters
n_pcs = 20 # Top dimensions to consider
top_freq_sites = 2E4 # Threshold to chose the most frequent sites across cells. Works well with 20,000
inputdir = 'cellcoverage'

if (!is.null(theargs$npcs)) {
   n_pcs <- as.integer(theargs$npcs)
}
if (!is.null(theargs$nclades)) {
   N_clades <- as.integer(theargs$nclades)
}
if (!is.null(theargs$ntfs)) {
   top_freq_sites <- as.integer(theargs$ntfs)
}
if (!is.null(theargs$inputdir)) {
   inputdir <- theargs$inputdir
}

print(sprintf("Using %d principal components", n_pcs))
print(sprintf("Picking %d clades", N_clades))
print(sprintf("Using %d top frequent sites", top_freq_sites))


#### Libs
library(data.table)
library(dplyr)
library(Matrix)
library(BuenColors)
library(stringr)
library(cowplot)
library(RColorBrewer)
library(ggpubr)
library(irlba)
library(umap)
library(gplots)
source('cusanovich2018_lib.r')


  # inputs
outdir = 'results/'



        # set directories
        system(paste0('mkdir -p ', outdir, '/data'))
        system(paste0('mkdir -p ', outdir, '/data/'))
        system(paste0('mkdir -p ', outdir, '/general'))
        pdf(paste0(outdir, 'general/', 'P3_identify_clades.pdf'))



#### Load window coverage per cell into a matrix
## Read genome windows and merge
df_regions = read.csv(file = 'genome_w5kbed', sep = '\t', header = F, stringsAsFactors = F) %>% unique(.)

winnames = paste(df_regions$V1, df_regions$V2, df_regions$V3, sep = '_')

print(head(winnames))

## List per cell covarage files
files = list.files(inputdir, pattern = '\\.txt$', full.names = F)
# Define cell barcode names
cell_barcodes = sapply(strsplit(files, '\\.'), "[", 1)
## Summary
message('Number of cells: ', length(files))
message('Number of regions: ', nrow(df_regions))



#### Build sparce empty matrix - This can be optimized!!!!
message('Load individual coverage files ...')
M = Matrix(
    0,
    nrow = length(winnames),
    ncol = length(cell_barcodes),
    sparse = T,
    dimnames = list(winnames, cell_barcodes))

## Fill values from coverage files
# Note we are only loading coverage for regions that are not V4=0 in the sample
for ( i in 1:length(files)) {
  f = paste0('cellcoverage/', files[i])
  message(f, ' ', i)
  d = fread(f, header = F, nThread = 4)
  idx = paste(d$V1, d$V2, d$V3, sep = '_')
  barcode = cell_barcodes[i]
  M[idx, barcode] = d$V4
}


# hierverder. writing matrix does not work.  actual goal is to insepct this error:
#  TF-IDF transformation
#  Dimension reduction
#  Error in irlba(tf_idf_counts, n_pcs, n_pcs) : 
#    starting vector near the null space
#  Execution halted

# the object looks like this:
#242506 x 100 Matrix of class "dgeMatrix"
#                              w5k w5k w5k w5k w5k w5k w5k w5k w5k w5k w5k w5k
#chr1_10000_15000     1.903095e-05 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
#chr1_15000_20000     1.903095e-05 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
#chr1_20000_25000     1.903095e-05 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
#chr1_25000_30000     1.903095e-05 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
#chr1_30000_35000     1.903095e-05 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
#chr1_35000_40000     1.903095e-05 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
#chr1_40000_45000     1.903095e-05 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
#chr1_55000_60000     1.903095e-05 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
# this is pbb the tf_idf_counts object.


# write.table(M, file=paste0(outdir, 'data/thetable.txt'), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
print("Stage 1 M loading")


# ## Binarize matrix
message('Binarize')
binary_mat = (M > 0) + 0
saveRDS(binary_mat, file = 'cus_P3_M.rds')
binary_mat       = readRDS('cus_P3_M.rds')
n_cells_with_site = rowSums(binary_mat)
options(repr.plot.width = 8, repr.plot.height = 4)
par(mfrow = c(1, 2))
hist(log10(n_cells_with_site), main = 'No. of Cells Each Site is Observed In', breaks = 50)


print("Stage 2 M loading")


# TODO: catch when there are not enough sites.
#### Filter bins
# Retain the most commonly used sites (top 50,000 here)
print(dim(binary_mat))
f_binary_mat =  binary_mat[ which(
        n_cells_with_site >= n_cells_with_site[order(n_cells_with_site, decreasing = T)[ top_freq_sites ]]
        ), ]
print(dim(f_binary_mat))
n_sites_per_cell = colSums(f_binary_mat)
f_binary_mat = f_binary_mat[rowSums(f_binary_mat) > 0, ]
print(dim(f_binary_mat))
print(sum(f_binary_mat > 0))
hist(log10(n_sites_per_cell), main = 'Number of Sites Each Cell Uses', breaks = 50)


print("Stage 3 Filter bins")


# #### Filter cells
# # Remove the lowest 10% of cells (in terms of site coverage) and ensure that there are now empty sites.
# f_binary_mat = f_binary_mat[, n_sites_per_cell >= 0]
# dim(f_binary_mat)
# f_binary_mat = f_binary_mat[rowSums(f_binary_mat) > 0, ]
# dim(f_binary_mat)



#### Transform the data using TF-IDF (term frequency–inverse document frequency)
message('TF-IDF transformation')
nfreqs = t(t(f_binary_mat) / Matrix::colSums(f_binary_mat))
idf = as(log(1 + ncol(f_binary_mat) / Matrix::rowSums(f_binary_mat)), 'sparseVector')
tf_idf_counts = as(Diagonal(x = as.vector(idf)), 'sparseMatrix') %*% nfreqs
print(dim(tf_idf_counts))

# elbow plot
# p_elbow = elbow_plot(tf_idf_counts, num_pcs = 200, title = 'PCA')
# p_elbow


print("Stage 4 TF-IDF")



#### Dimension reduction with truncated SVD (these are the two primary steps of LSI)
message('Dimension reduction')
# Here, we only retain components 2-n_pcs (component 1 is highly correlated with read depth) 
# This step can take a minute
set.seed(0) #For reproducibility
print(tf_idf_counts)
print("before")
SVD = irlba(tf_idf_counts, n_pcs, n_pcs)
print("after")
sk_diag = matrix(0, nrow = n_pcs, ncol = n_pcs)
diag(sk_diag) = SVD$d
# remove component 1 as suggested
sk_diag[1, 1] = 0
SVDumap_vd = t(sk_diag %*% t(SVD$v))
dim(SVDumap_vd)
# and truncate the distribution of LSI values at +/-1.5.
LSI_out = t(t(sk_diag %*% t(SVD$v)) %*% t(SVD$u))
LSI_out = t(scale(t(LSI_out)))
LSI_out[LSI_out > 1.5] = 1.5
LSI_out[LSI_out < -1.5] = -1.5
dim(LSI_out)
hist(LSI_out)


print("Stage 5 Dim reduction")


#### Find clusters
message('Find clusters')
# This step can take a minute too
hclust_cells = hclust(proxy::dist(t(sk_diag %*% t(SVD$v)), method = 'cosine'), method = 'ward.D2')
hclust_genes = hclust(proxy::dist(t(sk_diag %*% t(SVD$u)), method = 'cosine'), method = 'ward.D2')
color_pal = colorRampPalette(brewer.pal(8, 'Set2'))(N_clades)
hmcols = colorpanel(100, 'steelblue', 'white', 'tomato')
cells_tree_cut = cutree(hclust_cells, N_clades)
lsi_cells = cbind(colnames(f_binary_mat), cells_tree_cut)
## Save data
df_pseudobulk = lsi_cells
write.table(df_pseudobulk, 'cus_P3_clades.tsv', quote = F, sep = '\t', row.names = F, col.names = F)
## heatmap
message('Plotting clusters')
options(repr.plot.width = 4, repr.plot.height = 6)
heatmap.2(LSI_out, 
          col = hmcols, 
          ColSideColors = color_pal[as.factor(cells_tree_cut)], 
          #RowSideColors = color_pal[as.factor(genes_tree_cut)], 
          Rowv = as.dendrogram(hclust_genes), Colv = as.dendrogram(hclust_cells), 
          labRow = F, labCol = F, trace = 'none', scale = 'none', 
          useRaster = T)
# umap
df_umap_LSI <- run_umap(t(SVDumap_vd))
df_umap_LSI = as.data.frame(df_umap_LSI)
row.names(df_umap_LSI) = colnames(binary_mat)
labels2 = as.character(cells_tree_cut)
n_colors = length(unique(labels2))
colormap2 = colorRampPalette(brewer.pal(8, 'Set2'))(n_colors)
names(colormap2) = as.character(seq(1, n_colors))
p_LSI_cluster <- plot_umap(df_umap_LSI, labels = labels2, colormap = colormap2, title = 'Cusanovich2018')
p_LSI_cluster = p_LSI_cluster+labs(color = 'cluster')
p_LSI_cluster

print("Stage 4 Clusters")


dev.off()


