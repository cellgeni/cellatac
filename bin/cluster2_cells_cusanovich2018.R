#### Script to identify clusters of cells for a per cluster-based peak calling
#
# Adapted from pinellolab/scATAC-benchmarking Github repo
#   > Extra/10x_PBMC_5k/test_pseudobulk/Cusanovich2018_10xpbmc5k_idenfify_clades.ipynb 
# Based on http://atlas.gs.washington.edu/fly-atac/docs/
#


theargs <- R.utils::commandArgs(asValues=TRUE)

N_clades = 10 # Number of expected cell clusters
n_pcs = 20 # Top dimensions to consider

fn_matrix = "no-such-matrix"
fn_regions = "no-regions-given"
fn_cells = "no-cells-whatsoever"

if (!is.null(theargs$npcs)) {
   n_pcs <- as.integer(theargs$npcs)
}
if (!is.null(theargs$nclades)) {
   N_clades <- as.integer(theargs$nclades)
}

fn_matrix   <- theargs$matrix
fn_regions  <- theargs$regions
fn_cells <- theargs$cells

if (is.null(fn_matrix) || is.null(fn_regions) || is.null(fn_cells)) {
  print("Use --matrix --regions and --cells to specify matrix in market format")
  stopifnot(false)
}


print(sprintf("Using %d principal components", n_pcs))
print(sprintf("Picking %d clades", N_clades))


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

pdf('P3_identify_clades.pdf')

#### Load window coverage per cell into a matrix
## Read genome windows and merge


print("Stage 1 and 2 M loading, summary statistics")

f_binary_mat <- readMM(file = fn_matrix)
message('foo')
regions.names = read.delim(fn_regions, header = FALSE, stringsAsFactors = FALSE)
message('fol')
cells.names = read.delim(fn_cells, header = FALSE, stringsAsFactors = FALSE)
message('fok')
colnames(f_binary_mat) = cells.names$V1
message('foi')
rownames(f_binary_mat) = regions.names$V1
message('bar')
#print(colSums(f_binary_mat))
#print(rowSums(f_binary_mat))

# fixme get rid of hardcoded names win.stats and cell20000.stats
winstats <- read.table("matrix/win.stats", row.names=1, check.names=F, sep="\t", header=T, as.is=T, quote="", comment.char="")
message('bar')
options(repr.plot.width = 8, repr.plot.height = 4)
par(mfrow = c(1, 2))
hist(log10(winstats$degree), main = 'No. of Cells Each Site is Observed In', breaks = 50)

cellstats <- read.table("matrix/cell20000.stats", row.names=1, check.names=F, sep="\t", header=T, as.is=T, quote="", comment.char="")
hist(log10(cellstats$degree), main = 'Number of Sites Each Cell Uses', breaks = 50)



print("Stage 3 Filter bins")

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
SVD = irlba(tf_idf_counts, n_pcs, n_pcs)
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
row.names(df_umap_LSI) = colnames(f_binary_mat)
labels2 = as.character(cells_tree_cut)
n_colors = length(unique(labels2))
colormap2 = colorRampPalette(brewer.pal(8, 'Set2'))(n_colors)
names(colormap2) = as.character(seq(1, n_colors))
p_LSI_cluster <- plot_umap(df_umap_LSI, labels = labels2, colormap = colormap2, title = 'Cusanovich2018')
p_LSI_cluster = p_LSI_cluster+labs(color = 'cluster')
p_LSI_cluster

print("Stage 4 Clusters")


dev.off()


