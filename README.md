# cellatac

1. [Introduction](#introduction)
2. [Basic workflow](#basic-workflow)
3. [Cellatac functionality](#cellatac-functionality)
4. [Running cellatac](#running-cellatac)
    - [Cellatac needs](#cellatac-needs)
    - [Useful options](#useful-options)
    - [Example invocations](#example-invocations)
5. [Downstream analysis](#downstream-analysis)
6. [Outputs](#outputs)

## Introduction

Sanger Cellular Genetics ATAC-seq pipeline by Luz Garcia Alonso, Ni Huang and
Stijn van Dongen.

**cellatac** takes scATAC-seq aligned data (such as the fragments file from
Cell Ranger ATAC) and outputs a _count matrix of accessible chromatin peaks by
cell_ (i.e. analogous to the `filtered_peak_bc_matrix` from Cell Ranger ATAC).
The output matrix can then be used for dowstream analysis in Seurat, Scanpy,
cisTopic or any other tool.

Cell Ranger ATAC identifies the peaks by aggregating the signal of all the
barcodes in the sample. There are some papers reporting that this may be
unsuitable to detect peaks appearing in rare cell types/states. **cellatac**
uses [Cusanovich approach](https://www.sciencedirect.com/science/article/pii/S0092867418308559)
to increase the peak detection sensitivity by, first, identifying cell clusters
on a windows x cell rather than peaks per cell matrix, and then doing a peak
calling for each cluster. 

## Basic workflow

1. **Compute window coverage**. The genome is broken into 5kb windows and then
  each cell is scored for insertions in each window, generating a binary matrix
  (large and sparce) of windows by cells. Note that if multiple samples are
  provided, these are aggregated into a unique matrix.
2. **Cluster cells based on window coverage**. Matrix is filtered to retain
  only top 200k most commonly used windows. Using `Signac`, the binary matrix is
  normalized with Term Frequency-Inverse Document Frequency (TF-IDF) approach
  followed by a dimensionality reduction step using Singular Value Decomposition
  (SVD). The first LSI component is ignored as it often captures sequencing depth
  (technical variation) rather than biological variation. The 2-30 top remaining
  components are used to perform graph-based Louvain clustering (at a X
  resolution) and clusters are reported.
3. **Accessible chromatin peak calling per cluster**. Peaks are called
  separately on each cluster using `macs2`.
4. **Merge per-cluster peaks and generate peak by cell matrix.** Peaks from all
  clusters are merged into a master peak set (i.e. overlapping peaks are
  aggregated), and the corresponfding peak by cell matrix (indicating any reads
  occuring in each peak for each cell) is reported. Note that if multiple samples
  are provided, these are aggregated into a unique matrix. **This is the relevant
  matrix that you should use for clustering.**

## Cellatac functionality

* The clustering approach from the Cusanovich 2018 manuscript.
* Joint analysing of multiple 10x samples.
* A clustering step utilising Seurat.
* User-specified clustering.
* Peak/cell matrix based on merging per-cluster peaks.
* Peak/cell matrix per-cluster.


## Running cellatac

### Cellatac needs

* Singularity


### Useful options

```
--mermul true           merge multiplets using CR bam file
--mermul false          [default] use CR fragments.tsv.gz

--usecls __seurat__        [default] use Seurat/Signac approach resembling Cusanovich. It uses Louvain clustering instead.
--usecls __cusanovich__    use cusanovich-strict approach. It uses bi-clustering of cells and windows based on cosine distances using the ward algorithm.
--usecls <filename>        use custom clustering

--mergepeaks true       [default] merge cluster peaks, compute master cell/peak matrix
--perclusterpeaks false [default] computer per-cluster cell/peak matrix  
                            Note both can be set to true.

--cellbatchsize 500     [default] parallelisation bucket size (number of cells per bucket)
--nclades 10            [default] number of clusters to use (only applies to cusanovich-strict approach)
--sampleid <tag>        use <tag> in naming outputs. Not yet consistently applied
```


### Example invocations

This pipeline will need a singularity installation.  It supports two executing
platforms, *local* (simply execute on the machine you're currently on) and
*lsf*. To use the latter specify `-profile lsf`.


```
source=cellgeni/cellatac

manifest=/some/path/to/singlecell.csv
posbam=/some/path/to/possorted_bam.bam
fragments=/some/path/to/fragments.tsv.gz

cellbatchsize=400
nclades=10

nextflow run $source        \
  --cellcsv $manifest       \
  --fragments $fragments    \
  --cellbatchsize $cellbatchsize   \
  --posbam $posbam          \
  --outdir results          \
  --sampleid CR12345678     \
  -profile local            \
  --mermul true             \
  --usecls __seurat__       \
  --mergepeaks true         \
  -with-report reports/report.html \
  -resume -w work -ansi-log false \
  -config my.config
```

where `my.config` supplies singularity mount options and tells nextflow how
many CPUs it can utilise when using the local executor, e.g.

```
singularity {
  runOptions = '-B /some/path1 -B /another/path2'
  cacheDir = '/home/jovyan/singularity/'
}

executor {
    cpus   = 56
    memory = 600.GB
}
```

To run multiple samples:

```
nextflow run $source        \
  --muxfile mux.txt         \
  --cellbatchsize $cellbatchsize   \
  --outdir results          \
  -profile local            \
  --usecls __seurat__       \
  --mermul false            \
  --mergepeaks true         \
  -with-report reports/report.html \
  -resume -w work -ansi-log false \
  -c my.config
```

where `mux.txt` is a tab separated file that looks like this:

```
1   sampleX   /path/to/cellranger/output/for/sampleX
2   sampleY   /path/to/cellranger/output/for/sampleY
3   sampleZ   /path/to/cellranger/output/for/sampleZ
4   sampleU   /path/to/cellranger/output/for/sampleU
```

The first column will be used to make the barcodes in each sample unique across
the merged samples. As such it can be anything, but it is suggested to simply
use a range of integers starting at 1, or to use the last one or two
significant digits of the sample ID provided they are unique to each sample.

The cellranger output directories need not contain the full output. Currently
the pipeline expects these files:

```
fragments.tsv.gz  possorted_bam.bam singlecell.csv
```

When running multiple samples, the bam file is only used for its header. It is
possible to substitute the original bam file with the output of `samtools view
-H possorted_bam.bam`. This can be useful if it is necessary to copy the data
prior to running this pipeline; it is not necessary in this case to copy the
full position sorted bam file (they tend to be very large).  Currently it is
necessary that the substituted file has the same name `possorted_bam.bam`.

## Downstream analysis

The snippet below shows how to read in cellatac output as a Seurat object.

```
### Load scATAC binary matrix
# This is analogous to the gene expression count matrix used to analyze single-cell RNA-seq. 
# However, instead of genes, each row of the matrix represents a PEAK of the genome learned by cellatac. 
# The matrix is not binary, > 0 if there is any Tn5 cut site for each single barcode (i.e. cell) that map within each peak.
f_binary_mat <- readMM(file = paste0(cellatac_dir, 'peak_matrix/peaks_bc_matrix.mmtx.gz'))
regions.names = read.delim(paste0(cellatac_dir, 'peak_matrix/peaks.txt'), header = FALSE, stringsAsFactors = FALSE)
cells.names = read.delim(paste0(cellatac_dir, 'peak_matrix/bc.txt'), header = FALSE, stringsAsFactors = FALSE)
colnames(f_binary_mat) = cells.names$V1
rownames(f_binary_mat) = regions.names$V1

# Make binary
f_binary_mat@x[f_binary_mat@x > 0] <- 1

### Get some stats
# check distributions
message('Matrix size:\n', 'rows ', f_binary_mat@Dim[1], '\ncolumns ', f_binary_mat@Dim[2])
n_cells_with_site = rowSums(f_binary_mat)
options(repr.plot.width = 8, repr.plot.height = 4)
par(mfrow = c(1, 2))
hist(log10(n_cells_with_site), main = 'No. of Cells Each Site is Observed In', breaks = 50)
hist(n_cells_with_site, main = 'No. of Cells Each Site is Observed In', breaks = 50)

sites_per_cell = colSums(f_binary_mat)
options(repr.plot.width = 8, repr.plot.height = 4)
par(mfrow = c(1, 2))
hist(log10(sites_per_cell), main = 'No. of Sites Observed per Cell', breaks = 50)
hist(sites_per_cell, main = 'No. of Sites Observed per Cell', breaks = 50)

# compare coverage vs peak length 
pos = sapply(strsplit(rownames(f_binary_mat), split= ':'), tail , 1)
pos_len = sapply(strsplit(pos, split= '-'), function(x) as.numeric(x[2])-as.numeric(x[1]) )
par(mfrow = c(1, 1))
plot(pos_len, n_cells_with_site)
abline(v = f_binary_mat@Dim[2]*0.75)


# filter non-informative peaks: length < 2k bp or >75% frequency
f_binary_mat = f_binary_mat[ pos_len < 2000 & n_cells_with_site < f_binary_mat@Dim[2]*0.75, ]


# create CreateSeuratObject
chrom_assay <- CreateChromatinAssay(
  counts = f_binary_mat,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = NULL,
  min.cells = 5,
  min.features = 100
)

so <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
```

## Outputs

The most import outputs are described below. Cellatac creates a toplevel
output directory by default called 'results' (change with the `--outdir` option).

```
results/qc/seurat-clades.tsv        (cluster annotation)
results/qc/seurat.pdf               (cluster-annotated and sample-annotated UMAP plots)
results/peak_matrix/bc.txt                    (barcode(cell) labels)
results/peak_matrix/peaks.txt                 (peak labels)
results/peak_matrix/peaks_bc_matrix.mmtx.gz   (main output object)
results/peak_matrix/bc_peaks_matrix.mmtx.gz   (transpose of above)
results/cellmetadata/singlecell.tsv           (joined metadata)
results/cellmetadata/tagmap.txt               (links two-digit sampletag and samplename)
```

The list of all directories with a short description:

```
cellmetadata      (see above)
peak_matrix       (see above)
qc                (see above)
clus_peak_matrix  (per-cluster results, multiple bundles as in peak_matrix)
macs2             (per-cluster macs2 results)
peaks             (bed files with peak information)
win_matrix        (outputs relating to the windows selected for clustering, can be ignored)
```

The bed files in `peaks` are `allclusters_peaks_sorted.bed` and
`allclusters_masterlist_sps.bed`.  The first is the concatenated and
position-sorted list of all per-cluster peaks.  The second is the result of
merging these peaks so that overlapping and book-ended peaks are joined in a
single peak. The infix `sps` indicates it is sample-position-sorted, indicating
we use the chromosome order as found in the cellranger bam file.


