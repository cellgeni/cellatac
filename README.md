# cellatac
Sanger Cellular Genetics ATAC-seq pipeline by Luz Garcia Alonso, Ni Huang and Stijn van Dongen.

##Basic workflow

* Optionally merge multiplets.
* Compute window coverage.
* Cluster cells based on window coverage.
* Run macs2 for each cluster.
* Compute per-cluster and/or global cell/peak matrix.

## Cellatac currently has/implements/supports

* Optional merging of multiplets.
* The clustering approach from the Cusanovich 2018 manuscript.
* A clustering step utilising Seurat.
* User-specified clustering.
* Peak/cell matrix based on merging per-cluster peaks.
* Peak/cell matrix per-cluster.

