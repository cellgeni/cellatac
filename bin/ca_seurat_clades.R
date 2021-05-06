library(Signac)
library(Seurat)
library('Matrix')

theargs <- R.utils::commandArgs(asValues=TRUE)
have_multisample  <- is.null(theargs$"single-sample")
do_clip <- is.null(theargs$noclip)      # with contortion apologies.

### NOTE/check: make sure CR singlecell.csv table works as expected if we merge multiplets ourselves

# cell.names  filtered_cell.stats  filtered_window_bc_matrix.mmtx.gz  regions.names  win.stats

### Load cellATAC windows binary matrix
f_binary_mat <- as(readMM(file = 'mmtx/filtered_window_bc_matrix.mmtx.gz'), "dgCMatrix")
regions.names = read.delim('mmtx/regions.names', header = FALSE, stringsAsFactors = FALSE)
cell.names = read.delim('mmtx/cell.names', header = FALSE, stringsAsFactors = FALSE)
colnames(f_binary_mat) = cell.names$V1
rownames(f_binary_mat) = regions.names$V1

### Load cell calling from cellranger - Seurat likes this
metadata <- read.csv(file = 'singlecell.csv', header = TRUE, row.names=1, sep="\t")

### Creating a Seurat object using the windows/cell matrix
so <- CreateSeuratObject(
  counts = f_binary_mat,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata
)

### Normalization -  term frequency-inverse document frequency (TF-IDF) 
# is a two-step normalization procedure, 
# that both normalizes across cells to correct for differences in cellular sequencing depth, 
# and across peaks to give higher values to more rare peaks
so <- RunTFIDF(so)

### Feature selection: 
# The largely binary nature of scATAC-seq data makes it challenging to perform ‘variable’ feature selection, as we do for scRNA-seq. 
# Instead, we can choose to use only the top n% of features (peaks) for dimensional reduction, 
# or remove features present in less that n cells with the FindTopFeatures function. 
# Here, we will all features, though we note that we see very similar results 
# when using only a subset of features (try setting min.cutoff to ‘q75’ to use the top 25% all peaks), with faster runtimes. 
# Features used for dimensional reduction are automatically set as VariableFeatures for the Seurat object by this function.
so <- FindTopFeatures(so, min.cutoff = 'q0')

### Dimensional reduction: 
# We next run a singular value decomposition (SVD) on the TD-IDF normalized matrix, 
# using the features (peaks) selected above. 
# This returns a low-dimensional representation of the object.

clipvalue = 1.5 # If we want to clip to reproduce Cusanovich approach.
if (!do_clip) {
  clipvalue = NULL
}
so <- RunSVD(
  object = so,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi',
  scale.max = clipvalue,
  fastpath=FALSE
)

### Non-linear dimension reduction 
so <- RunUMAP(object = so, reduction = 'lsi', dims = 1:30)
so <- FindNeighbors(object = so, reduction = 'lsi', dims = 1:30)

### Clustering
so <- FindClusters(object = so, verbose = FALSE)

### Plot clustering and sample information

pdf('seurat.pdf')
DimPlot(object = so, label = TRUE) + NoLegend()
if (have_multisample) {
  so@meta.data$sample_id = sapply(strsplit(colnames(so), split = '-'), head, 1)
  DimPlot(object = so, label = TRUE, group.by="sample_id")
}
dev.off()

# Save file
saveRDS(so, file = 'seuratObject_win_matrix.rds')
clades = as.data.frame(Idents(so))
write.table(clades, file =  'seurat-clades.tsv', col.names = F, quote = F, sep = '\t')


