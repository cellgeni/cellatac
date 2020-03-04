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

### Cellatac needs

* Singularity


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

--cellbatchsize 500     [default] parallelisation bucket size (number of cells per bucket)
--nclades 10            [default] number of clusters to use

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

where `my.config` supplies singularity mount options and tells nextflow how many CPUs it can utilise
when using the local executor, e.g.

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

The first column will be used to make the barcodes in each sample unique across the merged samples. As
such it can be anything, but it is suggested to simply use a range of integers starting at 1.


