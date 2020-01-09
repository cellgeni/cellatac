

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
params.nbcest        =  10000           // estimated number of non-empty barcodes.
params.nclades       =  10
params.ntfs          =  20000
params.npcs          =  20

params.mermul        =  false
params.usecls        =  false


if (!params.fragments || !params.cellcsv || !params.posbam) {
  exit 1, "Please supply --fragments <CR-fragment-file> --cellcsv <CR-cellcsv-file> --posbam <CR-posbam-file>"
}


ch_fragments_cr = params.mermul ? Channel.empty() : Channel.fromPath(params.fragments)
ch_usercls = params.usecls ? Channel.empty() : Channel.fromPath(params.usecls)

thecellfile = file(params.cellcsv)
thebamfile  = file(params.posbam)


      // fixme noteme stick idxstats in here if needed
process prepare_cr {

  tag "cr-prep $cellbatchsize"

  publishDir "${params.outdir}", pattern: 'cellmetadata',   mode: 'copy'

  when: !params.mermul

  input:
  set file(f_cells), file(posbam) from Channel.from([[thecellfile, thebamfile]])
  val cellbatchsize from params.cellbatchsize
  val winsize from params.winsize

  output:
  file('cellmetadata/cells.tab')    into ch_celltab_cr
  file('cellmetadata/win.tab')      into ch_wintab_cr
  file('cellmetadata/cells.names')  into ch_index_names_cr
  file('cellmetadata')              into ch_metadata_cr
  file('cellmetadata/sample.chrlen')    into (ch_chrom_length_cr, ch_chrom_length2_cr)
  file('c_c.*') into ch_demux_batch_cr

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


process prepare_mm {        // merge multiplets

  tag "mm-prep $cellbatchsize"

  publishDir "${params.outdir}", pattern: 'cellmetadata_mm',   mode: 'copy'

  when: params.mermul

  input:
  set file(f_cells), file(posbam) from Channel.from([[thecellfile, thebamfile]])
  val cellbatchsize from params.cellbatchsize
  val winsize from params.winsize

  output:
  file('cellmetadata/cells.tab')    into ch_celltab_mm
  file('cellmetadata/win.tab')      into ch_wintab_mm
  file('cellmetadata/cells.names')  into ch_index_names_mm
  file('cellmetadata')              into ch_metadata_mm
  file('cellmetadata/sample.chrlen')    into (ch_chrom_length_mm, ch_chrom_length2_mm)
  file('c_c.*') into ch_demux_batch_mm

  file('fragmints2.tsv.gz')         into ch_fragments_mm

  shell:
  '''
  mkdir -p cellmetadata
# Chrosome length file.
  samtools view -H !{posbam}  \\
    | grep '@SQ'$'\\t''SN:'    \\
    | perl -ne '/\\bSN:(\\S+)/ && ($name=$1); /\\bLN:(\\d+)/ && ($len=$1); print "$name\\t$len\\n";' \\
    | uniq                    \\
    > cellmetadata/sample.chrlen.all

# get the main chromosomes. WARNING DANGERSIGN very crude regular expression filter.
# This filter basically avoids underscores and allows otherwise alphanumerical.
# TODO print joined string of all chromosomes for user perusal.
  grep -i 'chr[a-z0-9][a-z0-9]*\\>' cellmetadata/sample.chrlen.all > cellmetadata/sample.chrlen

#
  extract-fragments !{posbam} fragmints.tsv.gz !{task.cpus}

#
  merge-multiplets --chrom cellmetadata/sample.chrlen --max-n !{params.nbcest} --debug --outfrg fragmints2.tsv.gz fragmints.tsv.gz cellmetadata/bc-mm.map

# Just the names of selected cells.
  cut -f 2 cellmetadata/bc-mm.map | sort -u > cellmetadata/cells.names

# Tab file for mcx matrix loading
  nl -v0 -nln -w1 < cellmetadata/cells.names > cellmetadata/cells.tab

# Batch lists for demuxing
  split -l !{cellbatchsize} cellmetadata/cells.names c_c.

# Tab file for windows
  ca_make_chromtab.pl !{winsize} cellmetadata/sample.chrlen > cellmetadata/win.tab
  '''
}


// The merge-multiplet code was inserted in a clone of the preparation process.
// below are all the merging Y-junctions for plumbing the channels to the same destination.

  ch_celltab_cr.mix(ch_celltab_mm).set { ch_celltab }
  ch_wintab_cr.mix(ch_wintab_mm).set { ch_wintab }
  ch_index_names_cr.mix(ch_index_names_mm).set { ch_index_names }
  ch_metadata_cr.mix(ch_metadata_mm).set { ch_metadata }
  ch_chrom_length_cr.mix(ch_chrom_length_mm).set { ch_chrom_length }
  ch_chrom_length2_cr.mix(ch_chrom_length2_mm).set { ch_chrom_length2 }
  ch_demux_batch_cr.mix(ch_demux_batch_mm).set { ch_demux_batch }

  ch_fragments_cr.mix(ch_fragments_mm).set { ch_fragments }


/* The demux part does a complicated thing with buckets.
   This is a remnant from a time when everything was published in a bucket structure first.
   It will probably be removed.
*/

process demux {

  tag "$thetag"

  input:
  file cells from ch_demux_batch.flatten()
  file frags from ch_fragments

  output:
  file('celldata/[ACGT][ACGT][ACGT][ACGT]/*.bed') into ch_demuxed
  // file('*.list')  into ch_index
  file('celldata/mtx*.edges') into ch_matrix

  shell:
  thetag = cells.toString() - 'c_c.'
  '''
  mkdir -p celldata
  zcat !{frags} | samdemux.pl --barcodefile=!{cells} --outdir=celldata --bucket --fragments --ntest=0 --fnedges=mtx.!{thetag}.edges
  # find celldata -name '*.bed' > !{thetag}.list
  '''
}


ch_demuxed
  .flatMap()
  .map { it -> it.toString() - ~/.*celldata\/[ACGT]{4}\// - ~/\.bed/ + '\t' + it.toString()  }
  .into { ch_cellpaths_cluster; ch_cellpaths_peakcov }
  

process make_big_matrix {

  container 'quay.io/cellgeni/cellclusterer'

  publishDir "${params.outdir}/cellmetadata", mode: 'link', pattern: 'cell2win.mcx'

  input:
  file all_edges from ch_matrix.collect()
  file celltab from ch_celltab
  file wintab  from ch_wintab

  output:
  file 'cell2win.mcx' into ch_cell2win
  set file('cell2win.mcx'), file(celltab), file(wintab) into ch_dump_big_matrix

  shell:
  '''
  cat !{all_edges} | mcxload --stream-split -abc - -strict-tabc !{celltab} -strict-tabr !{wintab} --write-binary -o cell2win.mcx
  sleep 3
  '''
}


process mmtx_big_matrix  {

  tag "raw_window_bc_matrix"

  container 'quay.io/cellgeni/cellclusterer'

  publishDir "${params.outdir}/win_matrix", mode: 'link'

  output:
  file 'raw_window_bc_matrix.mmtx.gz'
  file 'raw_*.txt'

  input:
  set file('cell2win.mcx'), file('celltab'), file('wintab') from ch_dump_big_matrix

  script:
  '''
  # fixme build new image with datamash
  # n_entries=$(mcx query -imx cell2win.mcx | tail -n +2 | datamash sum 2
  n_entries=$(mcx query -imx cell2win.mcx | tail -n +2 | perl -ane 'chomp; $sum+=$_; END { print $sum }')
  ca_make_mmtx.sh -r wintab -c celltab -m cell2win.mcx -e $n_entries -t integer -o raw_window_bc_matrix.mmtx.gz
  cp celltab raw_bc.txt
  cp wintab  raw_window.txt
  '''
}


//      ^----- split between demux part and analysis part -----_      //


process filter_big_matrix {

  tag "bed-mcx-mmtx"

  container = 'quay.io/cellgeni/cellclusterer'

  publishDir "$params.outdir/win_matrix", mode: 'link', pattern: 'other_publish/filtered*',
    saveAs: { fname -> fname - ~/other_publish\// }

  when: true

  input:
  file('cellmetadata') from ch_metadata
  file('cell2win.mcx') from ch_cell2win
  val ntfs      from  params.ntfs

  output:
  file('outputs') into ch_load_mmtx
  file('other_publish/filtered*')

  shell:          
  '''
              # fixme, filtered_window_bc_matrix is still an implicit output from ca_top_region.sh
  mkdir matrix
  mkdir outputs
  mkdir other_publish
  cd matrix
  ca_top_region.sh \\
      -c ../cellmetadata/cells.tab    \\
      -w ../cellmetadata/win.tab      \\
      -m ../cell2win.mcx              \\
      -n !{ntfs}                      \\
      -C ../outputs/filtered_cell.stats          \\
      -W ../outputs/win.stats                    \\
      -R ../outputs/regions.names                \\
      -X ../outputs/filtered_window_bc_matrix.mmtx.gz  \\
      -N ../outputs/cells.names

  ln ../outputs/regions.names    ../other_publish/filtered_win.txt
  ln ../outputs/cells.names      ../other_publish/filtered_bc.txt
  ln ../outputs/filtered_window_bc_matrix.mmtx.gz   ../other_publish
  ln ../outputs/filtered_cell.stats                 ../other_publish
  '''
}


process do_the_clustering {

  tag "cusanovich2018"

  container = 'quay.io/cellgeni/cellclusterer'

  publishDir "$params.outdir/qc", mode: 'link', pattern: 'cus.*'

  when: !params.usecls

  input:
  val nclades   from  params.nclades
  val sampleid  from  params.sampleid
  val npcs      from  params.npcs
  file('inputs') from ch_load_mmtx

  output:
  file('cus.obj.*.clades.tsv') into ch_P4_clades
  file('cus.qc.*.pdf')
  file('cus.obj.*')

  shell:
  '''
  ln -s !{baseDir}/bin/cusanovich2018_lib.r .
  R --slave --quiet --no-save --args  \\
  --nclades=!{nclades}                \\
  --npcs=!{npcs}                      \\
  --matrix=inputs/filtered_window_bc_matrix.mmtx.gz  \\
  --regionnames=inputs/regions.names  \\
  --cellnames=inputs/cells.names      \\
  --winstats=inputs/win.stats         \\
  --cellstats=inputs/filtered_cell.stats  \\
  --sampleid=!{sampleid}              \\
  < !{baseDir}/bin/cluster_cells_cusanovich2018.R
  '''
}

ch_usercls.mix(ch_P4_clades).set{ ch_clustering }

process clusters_index {

  tag "${cladefile}"

  input:
  file metafile from ch_cellpaths_cluster.collectFile(name: 'cov.inputs', newLine: true)
  file cladefile from ch_clustering

  output:
  file('clusinfo.cl*') into ch_clusterbam

  shell:
  '''
  splitbyclus.pl !{cladefile} !{metafile} clusinfo.cl
  '''
}


  // This process derives the cluster tag from the file name. We only need to do this once,
  // in subsequent processes we can pass this tag on.


process clusters_makebam {

  tag "${clustag}"

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


process clusters_macs2 {

  tag "${clustag}"

  container = 'fooliu/macs2'

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
    .map { it.join('\n') + '\n' }

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

  container = 'quay.io/cellgeni/cellclusterer'

  publishDir "${params.outdir}/peak_matrix", mode: 'link'


  input:
  file metafile from ch_cellpeak.flatMap { ls -> ls.collect{ it.toString() } }.collectFile(name: 'peak.inputs', newLine: true)
  file('masterpeak.bed') from ch_masterbed_sps2.collect()
  file('cellnames.txt') from ch_index_names

  output:
    // file('cell2peak.gz')
  file('peaks_bc_matrix.mmtx.gz')
  file('bc_peaks_matrix.mmtx.gz')
  file('peaks.txt')
  file('bc.txt')

  shell:
  '''
                  # fixme: hardcoded output names.
  ca_peak_matrix.sh -c cellnames.txt -p masterpeak.bed -i peak.inputs
  '''
}





