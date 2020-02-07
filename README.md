# cellatac
Sanger Cellular Genetics ATAC-seq pipeline by Luz Garcia Alonso, Ni Huang and Stijn van Dongen.

### Basic workflow

* Optionally merge multiplets.
* Compute window coverage.
* Cluster cells based on window coverage.
* Run macs2 for each cluster.
* Compute per-cluster and/or global cell/peak matrix.

### Cellatac currently has/implements/supports

* Optional merging of multiplets.
* The clustering approach from the Cusanovich 2018 manuscript.
* A clustering step utilising Seurat.
* User-specified clustering.
* Peak/cell matrix based on merging per-cluster peaks.
* Peak/cell matrix per-cluster.

### Useful options

```
--mermul true           merge multiplets using CR bam file
--mermul false          [default] use CR fragments.tsv.gz

--usecls __seurat__        [default] custom approach resembling ...
--usecls __cusanovich__    use cusanovich approach
--usecls <filename>        use custom clustering

--mergepeaks true       [default] merge cluster peaks, compute master cell/peak matrix
--perclusterpeaks false [default] computer per-cluster cell/peak matrix  
                            Note both can be set to true.

--ncell 0               [default] set e.g. to 1000 to test-run smaller number of cells
--cellbatchsize 500     [default] parallelisation bucket size (number of cells per bucket)
--nclades 10            [default] number of clusters to use

--sampleid <tag>        use <tag> in naming outputs. Not yet consistently applied
```


### Example invocation

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
  --nclades $nclades        \
  --outdir results          \
  --sampleid CR12345678     \
  -profile local            \
  --mermul true             \
  --usecls __seurat__       \
  --mergepeaks true         \
  --perclusterpeaks true    \
  -with-report reports/report.html \
  -resume -w work -ansi-log false
```


