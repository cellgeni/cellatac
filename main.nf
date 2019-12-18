

params.fragments     =  null            // CR fragments file.
params.posbam        =  null            // CR possorted bam file.
params.ncell         =  0               // Can be set for testing purposes; suggest using a
                                        // sampled fragment file, e.g.
// zcat fragments.tsv.gz | perl -ne 'print if $. % 11 == 0;' | gzip > fragments11.tsv.gz

params.cellcsv       =  null            // CR singlecell.csv file
params.cellbatchsize =  500             // for demuxing and for masterpeak coverage

params.outdir        =  'results'
params.sampleid      =  'thesamp'       // Used in some names of output files

params.winsize       =  5000            // bin size for genome
params.nclades       =  10
params.ntfs          =  20000
params.npcs          =  20


if (!params.fragments || !params.cellcsv || !params.posbam) {
  exit 1, "Please supply --fragments <CR-fragment-file> --cellcsv <CR-cellcsv-file> --posbam <CR-posbam-file>"
}


// Channel.fromPath(params.cellfile).until { false}.set { ch_get_cells }
thecellfile = file(params.cellcsv)
thefragfile = file(params.fragments)
thebamfile  = file(params.posbam)


      // fixme noteme stick idxstats in here if needed
process prepare {

  tag "$cellbatchsize"

  cpus 1
  memory 1.GB

  publishDir "${params.outdir}", pattern: 'cellmetadata',   mode: 'copy'

  input:
  set file(f_cells), file(fragments), file(posbam) from Channel.from([[thecellfile, thefragfile, thebamfile]])
  val cellbatchsize from params.cellbatchsize
  val winsize from params.winsize

  output:
  file('cellmetadata/cells.tab')    into ch_celltab
  file('cellmetadata/win.tab')      into ch_wintab
  file('cellmetadata/cells.names')  into ch_index_names
  file('cellmetadata')              into ch_metadata
  file('cellmetadata/sample.chrlen')    into (ch_chrom_length, ch_chrom_length2)
  file('c_c.*') into ch_demux_batch
                                    // file('chosen_cells.info')

  shell:
  filter = params.ncell > 0 ? "head -n ${params.ncell}" : "cat"
  '''
  mkdir -p cellmetadata
# Chrosome length file.
  samtools view -H !{posbam}  \\
    | grep '@SQ'$'\\t''SN:'    \\
    | perl -ne '/\\bSN:(\\S+)/ && ($name=$1); /\\bLN:(\\d+)/ && ($len=$1); print "$name\\t$len\\n";' \\
    | uniq                    \\
    > cellmetadata/sample.chrlen

# Names + info of selected cells.
  perl -F, -ane 's/,/\t/g; print if $F[9] == 1' !{f_cells} \\
     | (sort -rnk 2 || true) | !{filter} > cellmetadata/chosen_cells.info

# Just the names of selected cells. 
  cut -f 1 cellmetadata/chosen_cells.info | sort > cellmetadata/cells.names

# Tab file for mcx matrix loading
  nl -v0 -nln -w1 < cellmetadata/cells.names > cellmetadata/cells.tab

# Batch lists for demuxing
  split -l !{cellbatchsize} cellmetadata/cells.names c_c.

# Tab file for windows
  ca_make_chromtab.pl !{winsize} cellmetadata/sample.chrlen > cellmetadata/win.tab
  '''
}



/* The demux part does a complicated thing with buckets.
   This is a remnant from a time when everything was published in a bucket structure first.
   It will probably be removed.
*/

process demux {

  cpus 1
  memory 4.GB

  tag "$thetag"

  input:
  file cells from ch_demux_batch.flatten()
  file frag from thefragfile

  output:
  file('celldata/[ACGT][ACGT][ACGT][ACGT]/*.bed') into ch_demuxed
  // file('*.list')  into ch_index
  file('celldata/mtx*.edges') into ch_matrix

  shell:
  thetag = cells.toString() - 'c_c.'
  '''
  mkdir -p celldata
  zcat !{frag} | samdemux.pl --barcodefile=!{cells} --outdir=celldata --bucket --fragments --ntest=0 --fnedges=mtx.!{thetag}.edges
  # find celldata -name '*.bed' > !{thetag}.list
  '''
}


ch_demuxed
  .flatMap()
  .map { it -> it.toString() - ~/.*celldata\/[ACGT]{4}\// - ~/\.bed/ + '\t' + it.toString()  }
  .into { ch_cellpaths_cluster; ch_cellpaths_peakcov }
  

process matrix {
  cpus 1
  memory 20.GB

  publishDir "${params.outdir}/cellmetadata", mode: 'link'

  input:
  file all_edges from ch_matrix.collect()
  file celltab from ch_celltab
  file wintab  from ch_wintab

  output:
  file 'cell2win.mcx' into ch_cell2win

  shell:
  '''
  cat !{all_edges} | mcxload --stream-split -abc - -strict-tabc !{celltab} -strict-tabr !{wintab} --write-binary -o cell2win.mcx
  sleep 3
  '''
}


//      ^----- split between demux part and analysis part -----_      //


process clusters_define_cusanovich2018_P3_C {

  tag "bottleneck"

  publishDir "$params.outdir/qc", mode: 'link'

  when: true

  input:
  file('cellmetadata') from ch_metadata
  file('cell2winmtx') from ch_cell2win
  val nclades   from  params.nclades
  val npcs      from  params.npcs
  val sampleid  from  params.sampleid
  val ntfs      from  params.ntfs

  output:
  file('cus.obj.*.clades.tsv') into ch_P4_clades
  file('cus.qc.*.pdf')
  file('cus.obj.*')

  shell:          
  '''
  mkdir matrix
  cd matrix
  ca_top_region.sh \\
      -c ../cellmetadata/cells.tab    \\
      -w ../cellmetadata/win.tab      \\
      -m ../cell2winmtx               \\
      -n !{ntfs}                      \\
      -C filtered_cell.stats          \\
      -W win.stats                    \\
      -R regions.names                \\
      -N cells.names

  cd ..

  ln -s !{baseDir}/bin/cusanovich2018_lib.r .
  R --slave --quiet --no-save --args  \\
  --nclades=!{nclades}                \\
  --npcs=!{npcs}                      \\
  --matrix=matrix/mtx.gz              \\
  --winstats=matrix/win.stats         \\
  --cellstats=matrix/filtered_cell.stats  \\
  --sampleid=!{sampleid}              \\
  --regionnames=matrix/regions.names  \\
  --cellnames=matrix/cells.names      \\
  < !{baseDir}/bin/cluster_cells_cusanovich2018.R
  '''
}


process clusters_index_P4 {

  tag "${cladefile}"

  input:
  file metafile from ch_cellpaths_cluster.collectFile(name: 'cov.inputs', newLine: true)
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
  publishDir "${params.outdir}/clusdef", mode: 'link'

  input:
  set val(clustag), file(clusmetafile) from ch_clusterbam.flatMap().map { [ it.baseName - 'clusinfo.cl', it ] }

  output:
  set val(clustag), file("cluster.${clustag}.b??") into ch_clustermacs

  shell:
  '''
  cat !{clusmetafile} | tr '\\n' '\\0' > clusmetafile0
  sort -m -k 1,1V -k 2,2n --files0-from=clusmetafile0 | perl -pe 'chomp;$_.="\\t+\\n";' > cluster.!{clustag}.bed
  '''
}


process clusters_macs2_P4 {

  tag "${clustag}"

  publishDir "${params.outdir}/macs2", mode: 'link'

  input:
  set val(clustag), file(clusregionfile) from ch_clustermacs

  output:
  set file('*.xls'), file('*.bed')
  file('*.narrowPeak') into ch_combine_clusterpeaks

  shell:
  '''
  # source /nfs/cellgeni/miniconda3/bin/activate py2

  outdir=macs2.!{clustag}.out
  mkdir -p $outdir
  macs2 callpeak -t !{clusregionfile} -g 2.7e9 -n !{clustag} --outdir . -f AUTO \\
     --nomodel     \\
     --shift -100  \\
     --extsize 200 \\
     2> macs2.!{clustag}.log
  '''
}


process peaks_masterlist {

  tag "masterlist"

  publishDir "${params.outdir}/peaks", mode: 'link'

  input:
  file np_files from ch_combine_clusterpeaks.collect()
  file sample_idxstats from ch_chrom_length

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

  // publishDir "${params.outdir}/mp_counts"
  // A lot of files. This is simply the raw input for the cell/peak matrix.

  input:
  file(masterbed_sps) from ch_masterbed_sps.collect()
  file sample_chrlen from ch_chrom_length2.collect()

  file(celldef_list) from ch_cellpaths_peakcov
    .collate(params.cellbatchsize)
    .map { it.join('\n') }

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

  publishDir "${params.outdir}/peak_matrix", mode: 'link'

  input:
  file metafile from ch_cellpeak.flatMap { ls -> ls.collect{ it.toString() } }.collectFile(name: 'peak.inputs', newLine: true)
  file('masterpeak.bed') from ch_masterbed_sps2.collect()
  file('cellnames.txt') from ch_index_names

  output:
  file('cell2peak.gz')
  file('peaks.txt')
  file('cells.txt')

  shell:
  '''
  ca_peak_matrix.sh -c cellnames.txt -w masterpeak.bed -i peak.inputs
  '''
}



