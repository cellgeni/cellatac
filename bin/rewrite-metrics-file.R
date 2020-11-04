

theargs <- commandArgs(trailingOnly=TRUE)
fninput <- theargs[1]

if (is.null(fninput)) {
  print("supply <metrics-file> please (10x multiomis per_barcode_metrics.csv file)")
  print("Usage e.g. R --no-save --args e18_mouse_brain_fresh_5k_per_barcode_metrics.csv < rewrite-metrics-file.R")
  stopifnot(FALSE)
}

a <- read.table(fninput, check.names=F, sep="\t", header=T, as.is=T, quote="", comment.char="")

a$ATAC_passed_filters <- a$atac_raw_reads - a$atac_dup_reads - a$atac_chimeric_reads - a$atac_unmapped_reads -a$atac_mitochondrial_reads
a$ATAC_cell_id <- 'None'
a$ATAC_DNase_sensitive_region_fragments <- 2
a$ATAC_enhancer_region_fragments <- 2
a$ATAC_promoter_region_fragments <- 2
a$ATAC_on_target_fragments <- 2
a$ATAC_blacklist_region_fragments <- 2

cellcounter <- 0

for (i in seq_along(a$ATAC_cell_id)) {
  if (a$is_cell[i] > 0) {
    a$ATAC_cell_id[i] <- sprintf("_cell_%d", cellcounter)
    cellcounter <- cellcounter+1
  }
}

mymap <- list(
"barcode"="barcode",
"atac_raw_reads"="total",
"atac_dup_reads"="duplicate",
"atac_chimeric_reads"="chimeric",
"atac_unmapped_reads"="unmapped",
"atac_lowmapq"="lowmapq",
"atac_mitochondrial_reads"="mitochondrial",
"ATAC_passed_filters"="passed_filters",
"ATAC_cell_id"="cell_id",
"is_cell"="is__cell_barcode",
"atac_TSS_fragments"="TSS_fragments",

"ATAC_DNase_sensitive_region_fragments"="DNase_sensitive_region_fragments",
"ATAC_enhancer_region_fragments"="enhancer_region_fragments",
"ATAC_promoter_region_fragments"="promoter_region_fragments",
"ATAC_on_target_fragments"="on_target_fragments",
"ATAC_blacklist_region_fragments"="blacklist_region_fragments",

"atac_peak_region_fragments"="peak_region_fragments",
"atac_peak_region_cutsites"="peak_region_cutsites")

newnames <- names(mymap)

a2 <- a[,newnames]
colnames(a2) <- mymap[newnames]

write.table(a2, "singlecell2.csv", quote = F, sep = ',', row.names = F, col.names = T)

