

/* Sanger Cellular Genetics ATAC-seq pipeline.
 * Scripts developed by Luz Garcia Alonso
 * The first version impmlements the clustering approach from the Cusanovich 2018 manuscript.
 * Nextflow integration by Stijn van Dongen.
*/

/* Example input creation (barcode file creation omitted):
   samtools view -s 0.05 atac_v1_pbmc_10k_possorted_bam.bam > sample_pb.bam
   samtools view -h -o - sample_pb.bam | samdemux.pl --barcodefile=barcode100.txt --outdir=d100
   bedtools bamtobed -i sample_pb.bam  | cut -f 1-3 | uniq > sample_pb.bed
   cd d100
   parallel -j 8 sam2bam.sh ::: *-1.sam
   rm -f *.sam

     Then:
   full path name to  |   is argument to option
----------------------+--------------------------------------
   sample_pb.bam      |   --psbam
   sample_pb.bed      |   --psbed
   barcode100.txt     |   --cellfile
   d100               |   --cellbamdir
*/

/* NOTES/TODO/QUESTIONS

   - make samdemux.nf for demultiplexing, converting to bed, and creating fai

   - do we need f_psbed? See todo near it.

   - We now pass in the possorted bed file.
     For the psbed and psbam files, do we need to subset only those reads that belong to is__cell_barcode cells?

   ! these look nearly identical (check prepare step):
     *  possorted_bam.genome.txt
     *  possorted_bam.chromosomes.txt
          there is also the file
     *  hg38.chrom.size in assets; check how it's used, especially relative to sorting requirement.

     In the P6 scripts these seem to be made again, same as in P2 script.

   - the pipeline chunks cells into batches in a few places, using meta files
     of bam locations.
     The cluster R script stages all w5k.txt files in single directory;
     this can be changed to a meta file (constructed with collectFile).
     Need to change the R script input reading.

   - encode sample name in R script output.

   # Note: singularity used for R clustering process and for macs2.

   ? input tag to encode parameters, use in output dir? User could this outside pl.
   - chrEBV hardcoded
   ? improve inclusion of bin/cusanovich2018_lib.r (currently linked in)
   - more arguments for parameters?
   - make sure parameters and pipeline caching/resumption play as they should.
   - use caching for f_psbam outside NF perhaps, time consuming. StoreDir?
*/


params.chromsizes    =  "$baseDir/assets/hg38.chrom.size"
params.genome        =  'hg38'
params.outdir        =  'results'
params.cellbamdir    =   null
params.psbam         =   null
params.psbed         =   null
params.sampleid      =  'thesamp'
params.cellfile      =  null

params.cellbatchsize = 100            // some things parallelise over cells, but per-cell is overkill.

params.nclades       =  10
params.ntfs          =  2E4
params.npcs          =  20
params.windowsize    =  5000


if (!params.psbam || !params.psbed || !params.cellfile || !params.cellbamdir) {
  exit 1, "Please supply --psbam --psbed --cellfile --cellbamdir each with argument"
}

Channel.fromPath(params.cellfile).until { true}. set { ch_get_cells }

Channel.fromPath(params.cellfile).until { false}.set { ch_get_cells2 }


process genome_make_windows {

  tag  "${f_chromsizes}"

  publishDir "${params.outdir}/genome", method: 'copy'
  input:
  set val(gntag), file(f_chromsizes) from Channel.from([[params.genome, file(params.chromsizes)]])

  output:
  file '*.bed' into ch_genome_w5k, ch_genomebed_P3
  
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
  set val(gntag), val(sample), file(f_psbam), file(f_psbed) from Channel.from([[params.genome, params.sampleid, file(params.psbam), file(params.psbed)]])
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
  # 1b.2
  bedtools intersect -a !{gw5k} -b !{f_psbed} -wa | uniq > !{sample}.w5k.bed

  # 1b.3                                                                            #   !{sample}.idxstats
  samtools idxstats !{f_psbam} | cut -f 1-2 | uniq > !{sample}.idxstats

  # Order the windows bed file as it is in the chromosomes file
  # 1b.4                                                                            #   !{sample}.w5ksorted.bed
  bedtools sort -faidx !{sample}.idxstats -i !{sample}.w5k.bed | uniq > !{sample}.w5ksorted.bed
  '''
}

process cells_get_files {

    tag "${file_cellnames}"

    errorStrategy 'terminate'

    input:
        file file_cellnames from ch_get_cells2
    output:
        file 'themetafile' into ch_cellbams_chunk, ch_clusbams2

    shell:
    '''
    while read cellname; do
      bam="!{params.cellbamdir}/$cellname.bam"
      if [[ ! -e $bam ]]; then
        >&2 echo "Bam file $bam not found"
        # false
      else
        echo -e "$cellname\t$bam"
      fi
    done < !{file_cellnames} > themetafile
    '''
}

ch_cellbams_chunk
   .splitText(by: params.cellbatchsize, file: "celllist.")
   .into{ch_cellbams_coverage2; ch_cellbams_peaks2}


process cells_window_coverage_P2 {         /* P2 process cusanovich2018.P2.windowcount.sh */

  tag "${cellbam_list}"

  publishDir "$params.outdir/cells"

  input:
      file sample_w5ksorted from ch_sample_w5ksorted.collect()
      file sample_chrlen    from ch_chrom_length.collect()
      file cellbam_list from ch_cellbams_coverage2

  output:
    file('*.txt') into ch_cellcoverage_P3

  shell:
  '''
  while read cellname cellbam; do
    bedtools coverage -sorted -header \\
      -a !{sample_w5ksorted}      \\
      -b $cellbam                 \\
      -g !{sample_chrlen} | awk -F"\\t" '{if($4>0) print $0}' > $cellname.w5k.txt
  done < !{cellbam_list}
  '''
}


process clusters_define_cusanovich2018_P3 {

  tag "bottleneck"

  publishDir "$params.outdir/qc", pattern: '*.pdf', method: 'link'

  input:
  file('genome_w5kbed') from ch_genomebed_P3
  file('cellcoverage/*') from ch_cellcoverage_P3.flatMap().collect()

        // NOTE / TODO. we could generate a meta file and pass this to the R script.

  output:
  file('cus_P3_M.rds') into ch_Px_rds
  file('cus_P3_clades.tsv') into ch_P4_clades
  file('*.pdf')

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
  file metafile from ch_clusbams2
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
  publishDir "${params.outdir}/clusbam"

  input:
  set val(clustag), file(clusmetafile) from ch_clusterbam.flatMap().map { [ it.baseName - 'clusinfo.cl', it ] }

  output:
  set val(clustag), file('cluster*bam') into ch_clustermacs

  shell:
  '''
  samtools merge "cluster.!{clustag}.bam" -b !{clusmetafile}
  '''
}


process clusters_macs2_P4 {

  tag "${clustag}"

  publishDir "${params.outdir}/macs2"

  input:
  set val(clustag), file(clusbamfile) from ch_clustermacs

  output:
  set file('*.xls'), file('*.bed')
  file('*.narrowPeak') into ch_combine_clusterpeaks

  shell:
  '''
  # source /nfs/cellgeni/miniconda3/bin/activate py2

  outdir=macs2.!{clustag}.out
  mkdir $outdir
  macs2 callpeak -t !{clusbamfile} -f BAM -g 2.7e9 -n !{clustag} --outdir . \\
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
  file('allclusters_masterlist_sps.bed') into ch_masterbed_sps

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


process masterlist_cells_count {

  tag "${cellbam_list}"

  publishDir "${params.outdir}/mp_counts"

  input:
  file(masterbed_sps) from ch_masterbed_sps.collect()
  file sample_chrlen from ch_chrom_length2.collect()

  file(cellbam_list) from ch_cellbams_peaks2

  output:
  file('*.mp.txt')

  shell:
  '''
  while read cellname cellbam; do
    bedtools coverage               \\
    -a !{masterbed_sps}             \\
    -b $cellbam -sorted -header     \\
    -g !{sample_chrlen} | awk -F"\t" '{if($4>0) print $0}' > $cellname.mp.txt
  done < !{cellbam_list}
  '''
}



