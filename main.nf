

/* Sanger Cellular Genetics single cell ATAC-seq pipeline.
 * Scripts developed by Luz Garcia Alonso.
 * The first version impmlements the clustering approach from the Cusanovich 2018 manuscript.
 * Nextflow integration by Stijn van Dongen.
*/

/* NOTES/TODO/QUESTIONS

   ! cellnames duplicates first column of cellpaths. simplify.

   ! sort -k 1,1V assumes 'V' (human alphanumeric) is the input sort order (bam/fragment chromosome order);
     so it is assume that the fragment file and resulting per-cell fragment files are in that order,
     but we don't know for sure. It is not an unreasonable assumption and bedtools will probably
     complain if there is a mismatch somewhere.

     We use this assumpion when merging bed files to get a per-cluster bed file to feed to macs2.

   - cellatac scripts have implicit naming convention. fix.
   - chrEBV hardcoded
   - encode sample name in R script output.

   ! these look nearly identical (check prepare step):
     *  possorted_bam.genome.txt
     *  possorted_bam.chromosomes.txt
          there is also the file
     *  hg38.chrom.size in assets; check how it's used, especially relative to sorting requirement.

   # Note: singularity used for R clustering process and for macs2.

   ? input tag to encode parameters, use in output dir? User could do this outside pl.
   ? improve inclusion of bin/cusanovich2018_lib.r (currently linked in)
   - more arguments for parameters?
   - make sure parameters and pipeline caching/resumption play as they should.
   - use caching for f_psbam outside NF perhaps, time consuming. StoreDir?
*/


params.chromsizes    =  "$baseDir/assets/hg38.chrom.size"
params.genome        =  'hg38'
params.outdir        =  'results'

params.psbam         =   null
params.psbai         =   null
params.sampleid      =  'thesamp'
params.cellfile      =  null
params.cellnames     =  null

params.cellbatchsize = 100            // some things parallelise over cells, but per-cell is overkill.

params.nclades       =  10
params.ntfs          =  20000
params.npcs          =  20
params.windowsize    =  5000

params.fragments     = true          // aim to support both bam and bed/fragment input. Need to about it in
                                     // places. The cell file has format
                                     // <barcode|tag>  <path/to/file/that/is/bam/or/bed/demultiplexed/component>

NWIN = params.ntfs


if (!params.psbam || !params.cellfile) {
  exit 1, "Please supply --psbam --cellfile each with argument"
}

// Channel.fromPath(params.cellfile).until { false}.set { ch_get_cells }

thecellnames = file(params.cellnames)
thecellfile  = file(params.cellfile)


process genome_make_windows {

  tag  "${f_chromsizes}"

  publishDir "${params.outdir}/genome", method: 'copy'
  input:
  set val(gntag), file(f_chromsizes) from Channel.from([[params.genome, file(params.chromsizes)]])

  output:
  file '*.bed' into ch_genome_w5k, ch_genomebed_P3, ch_genomebed_P3_B
  
  shell:
  '''
  cat !{f_chromsizes} <(echo -e "chrEBV\t171823") > t.size
  bedtools makewindows -g t.size -w !{params.windowsize} >  !{gntag}.w5k.bed
  '''
}


process posbam_prepare_info {

  tag  "${f_psbam}"

  publishDir "$params.outdir/sample"

  input:
  set val(gntag), val(sample), file(f_psbam), file(f_psbai) from Channel.from([[params.genome, params.sampleid, file(params.psbam), file(params.psbai)]])
  file(gw5k) from ch_genome_w5k

  output:
  file "${sample}.w5ksorted.bed" into ch_sample_w5ksorted
  file "${sample}.chrlen" into (ch_chrom_length, ch_chrom_length2)   // called possorted_bam.chromosomes.txt in P2 script
  file "${sample}.idxstats" into ch_idxstats

  shell:
  '''
  # 1a
  samtools view -H !{f_psbam}   \\
    | grep '@SQ'$'\t''SN:'    \\
    | perl -ne '/\\bSN:(\\S+)/ && ($name=$1); /\\bLN:(\\d+)/ && ($len=$1); print "$name\\t$len\\n";' \\
    | uniq                    \\
    > !{sample}.chrlen                                                              #   !{sample}.chrlen
                                                                                    #   called possorted_bam.chromosomes.txt
                                                                                    #   in P2 script.

  # nf-NOTES perhaps useful to use a caching mechanism (outside NF), as this is very time-consuming
  # 1b.1                                                                            #   !{sample}.bed
  # bedtools bamtobed -i !{f_psbam} \\
  #   | cut -f 1-3 | uniq > !{sample}.bed
  # CHANGE: we now assume this has been done beforehand.

  # Intersect genome-wide windows to filter out regions that are not in the read data
  # TODO do we gain much from this? (measure how quick it is and how many windows we lose).
  # Now doing this. Original code first:
  # 1b.2
  # bedtools intersect -a !gw5k -b !f_psbed -wa | uniq > !sample.w5k.bed
  # ^^^^^^^^^^^^^^^^^^^^^^ original code, seems to take a long time.
  # See below for new code.
  # The old code led to a small reduction in number of windows (<10%).
  # The new code is much much faster. Empty windows should still be lost at the clustering stage.

  # If we simply copy the source to the destination:
  # cp !{gw5k} !{sample}.w5k.bed
  # Then there may be chromosome names in the first not occurring in the latter.
  # So, let's grep the chromosome names from the idxstats, which we'll compute first.

  # 1b.3                                                                            #   !{sample}.idxstats
  samtools idxstats !{f_psbam} | cut -f 1-2 | uniq > !{sample}.idxstats

  # 1b.2                                                                            #   !{sample}.idxstats
  grep -Fwf <(cut -f 1 !{sample}.idxstats) !{gw5k} > !{sample}.w5k.bed

  # Order the windows bed file as it is in the chromosomes file
  # 1b.4                                                                            #   !{sample}.w5ksorted.bed
  bedtools sort -faidx !{sample}.idxstats -i !{sample}.w5k.bed | uniq > !{sample}.w5ksorted.bed
  '''
}

// process cells_get_files {
// 
//     tag "${file_cellnames}"
// 
//     errorStrategy 'terminate'
// 
//     input:
//         file file_cellnames from thecellnames
//     output:
//         file 'themetafile' into ch_cellbams_chunk, ch_clusbams2
// 
//     shell:
//     '''
//     while read cellname; do
//       bam="!{params.cellbamdir}/$cellname.bam"
//       if [[ ! -e $bam ]]; then
//         >&2 echo "Bam file $bam not found"
//         # false
//       else
//         echo -e "$cellname\t$bam"
//       fi
//     done < !{file_cellnames} > themetafile
//     '''
// }

Channel.fromPath(params.cellfile)
   .splitText(by: params.cellbatchsize, file: "celllist.")
      .into{ch_celldefs_coverage2; ch_celldefs_peaks2}


process cells_window_coverage_P2 {         /* P2 process cusanovich2018.P2.windowcount.sh */

  tag "${celldef_list}"

  publishDir "$params.outdir/cells"

  input:
      file sample_w5ksorted from ch_sample_w5ksorted.collect()
      file sample_chrlen    from ch_chrom_length.collect()
      file celldef_list from ch_celldefs_coverage2

  output:
    file('*.txt') into (ch_cellcoverage_P3, ch_cellcoverage_P3_B)

  shell:
  '''
  while read cellname celldef; do
    bedtools coverage -sorted -header \\
      -a !{sample_w5ksorted}      \\
      -b $celldef                 \\
      -g !{sample_chrlen} | awk -F"\\t" '{if($4>0) print $0}' > $cellname.w5k.txt
  done < !{celldef_list}
  '''
}


process clusters_define_cusanovich2018_P3_B {

  tag "bottleneck"

  publishDir "$params.outdir/qc", pattern: '*.pdf', method: 'link'

  input:
  file('genome_w5kbed') from ch_genomebed_P3_B
  file('cellnames.txt') from thecellnames
  file metafile from ch_cellcoverage_P3_B.flatMap { ls -> ls.collect{ it.toString() } }.collectFile(name: 'cov.inputs', newLine: true)
  // file('cellcoverage/*') from ch_cellcoverage_P3_B.flatMap().collect()

  output:
  file('cus_P3_clades.tsv') into ch_P4_clades
  file('*.pdf')

  shell:          
  '''
  mkdir matrix
  cd matrix
  ca_top_region.sh -c ../cellnames.txt -w ../genome_w5kbed -i ../cov.inputs -n !{NWIN}
  # fixme: define outputs of ^ script using options. Currently implicit.
  cd ..

  ln -s !{baseDir}/bin/cusanovich2018_lib.r .
  R --slave --quiet --no-save --args  \\
  --nclades=!{params.nclades}         \\
  --npcs=!{params.npcs}               \\
  --matrix=matrix/mtx.gz              \\
  --regions=matrix/regions!{NWIN}.txt \\
  --cells=matrix/cells.txt            \\
  < !{baseDir}/bin/cluster2_cells_cusanovich2018.R
  '''
}


process clusters_define_cusanovich2018_P3 {

  tag "bottleneck"

  publishDir "$params.outdir/qc", pattern: '*.pdf', method: 'link'

  input:
  file('genome_w5kbed') from ch_genomebed_P3
  file('cellcoverage/*') from ch_cellcoverage_P3.flatMap().collect()

//  This script is now a dead-end; we retain it for now to make sure new approach does
//  the same.
//  output:
//  file('cus_P3_clades.tsv') into ch_P4_clades
//  file('*.pdf')

  shell:
  '''
  ln -s !{baseDir}/bin/cusanovich2018_lib.r .
  R --slave --quiet --no-save --args  \\
  --nclades=!{params.nclades}         \\
  --npcs=!{params.npcs}               \\
  --ntfs=!{params.ntfs}               \\
  --inputdir=cellcoverage             \\
  < !{baseDir}/bin/cluster_cells_cusanovich2018.R
  '''
}


process clusters_index_P4 {

  tag "${cladefile}"
            // fixme using old idiom for channels, but now have a list. as a hack I flatmap it.
            // so that I can then use collectFile.
  input:
  file metafile from thecellfile
  file cladefile from ch_P4_clades

  output:
  file('clusinfo.cl*') into ch_clusterbam

  shell:
  '''
  splitbyclus.pl !{cladefile} !{metafile} clusinfo.cl
  '''
}


  // This process derives the cluster tag from the file name. We only need to do this once,
  // in subsequent processes we can pass this tag on.


process clusters_makebam_P4 {

  tag "${clustag}"

  /* TODO: insert clade,pcs,ntfs in output directory name? */
  publishDir "${params.outdir}/clusdef"

  input:
  set val(clustag), file(clusmetafile) from ch_clusterbam.flatMap().map { [ it.baseName - 'clusinfo.cl', it ] }

  output:
  set val(clustag), file("cluster.${clustag}.b??") into ch_clustermacs

  shell:
  if (params.fragments) {
                    // sort 1,1V sorts alphanumerically, so chr9 chr10 rather than chr1 chr10.
                    // There is an assumption here that this is the input sort order.
    '''
    cat !{clusmetafile} | tr '\\n' '\\0' > clusmetafile0
    sort -m -k 1,1V -k 2,2n --files0-from=clusmetafile0 | perl -pe 'chomp;$_.="\\t+\\n";' > cluster.!{clustag}.bed
    '''
  }
  else {
    '''
    samtools merge "cluster.!{clustag}.bam" -b !{clusmetafile}
    '''
  }
}


process clusters_macs2_P4 {

  tag "${clustag}"

  publishDir "${params.outdir}/macs2"

  input:
  set val(clustag), file(clusregionfile) from ch_clustermacs

  output:
  set file('*.xls'), file('*.bed')
  file('*.narrowPeak') into ch_combine_clusterpeaks

  shell:
  fmtoption = params.fragments ? "-f AUTO" : "-f BAMPE"
  '''
  # source /nfs/cellgeni/miniconda3/bin/activate py2

  outdir=macs2.!{clustag}.out
  mkdir -p $outdir
  macs2 callpeak -t !{clusregionfile} -g 2.7e9 -n !{clustag} --outdir . !{fmtoption}\\
     --nomodel     \\
     --shift -100  \\
     --extsize 200 \\
     2> macs2.!{clustag}.log
  '''
}


process peaks_masterlist {

  tag "masterlist"

  publishDir "${params.outdir}/peaks"

  input:
  file np_files from ch_combine_clusterpeaks.collect()
  file sample_idxstats from ch_idxstats

  output:
  file('allclusters_peaks_sorted.bed')
  file('allclusters_masterlist_sps.bed') into (ch_masterbed_sps, ch_masterbed_sps2)

    // NOTE may want to encode some cluster parameters in the file name? Also preceding processes
  shell:
  '''
  cat !{np_files} | cut -f 1-3 | sort -k1,1 -k2,2n > allclusters_peaks_sorted.bed
  # use -d -1 to avoid mergeing regions overlapping only 1bp
  bedtools merge -i allclusters_peaks_sorted.bed -d -1 > allclusters_masterlist.bed

  # Note moved this step from P6 to here. sps == sample-pos-sorted, sorted according to sample bam.
  bedtools sort -faidx !{sample_idxstats} -i allclusters_masterlist.bed > allclusters_masterlist_sps.bed
  '''
}


process cells_masterlist_coverage {

  tag "${celldef_list}"

  publishDir "${params.outdir}/mp_counts"

  input:
  file(masterbed_sps) from ch_masterbed_sps.collect()
  file sample_chrlen from ch_chrom_length2.collect()

  file(celldef_list) from ch_celldefs_peaks2

  output:
  file('*.mp.txt') into ch_cellpeak

  shell:
  '''
  while read cellname celldef; do
    bedtools coverage               \\
    -a !{masterbed_sps}             \\
    -b $celldef -sorted -header     \\
    -g !{sample_chrlen} | awk -F"\t" '{if($4>0) print $0}' > $cellname.mp.txt
  done < !{celldef_list}
  '''
}



process make_peakmatrix {

  tag "peak-cell-matrix"

  publishDir "${params.outdir}/peak_matrix"

  input:
  file metafile from ch_cellpeak.flatMap { ls -> ls.collect{ it.toString() } }.collectFile(name: 'peak.inputs', newLine: true)
  file('masterpeak.bed') from ch_masterbed_sps2.collect()
  file('cellnames.txt') from thecellnames

  output:
  file('cell2peak.gz')
  file('peaks.txt')
  file('cells.txt')

  shell:
  '''
  ca_peak_matrix.sh -c cellnames.txt -w masterpeak.bed -i peak.inputs
  '''
}



