

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
params.ntfs_ofs      =  1
params.npcs          =  20
params.sumrank		   =  'true'

params.mermul        =  false
params.usecls        =  '__seurat__'
params.mergepeaks    =  true
params.perclusterpeaks  =  false
params.stopatcluster =  false

params.muxfile       =  null

if ((!params.fragments || !params.cellcsv || !params.posbam) && !params.muxfile) {
  exit 1, "Please supply --fragments <CR-fragment-file> --cellcsv <CR-cellcsv-file> --posbam <CR-posbam-file>"
}


ch_usercls = params.usecls =~ /^__.*__$/ ? Channel.empty() : Channel.fromPath(params.usecls)
ch_mux     = params.muxfile && !params.mermul ? Channel.fromPath(params.muxfile) : Channel.empty()


ch_cellfile = params.cellcsv    ? Channel.fromPath(params.cellcsv)   : Channel.empty()
ch_fragfile_cr = params.fragments  ? Channel.fromPath(params.fragments) : Channel.empty()
ch_bamfile  = params.posbam     ? Channel.fromPath(params.posbam)    : Channel.empty()

ch_cellfile.into { ch_cellfile_mm; ch_cellfile_cr; ch_cellfile_seurat }
ch_bamfile.into  { ch_bamfile_mm; ch_bamfile_cr }

/*
thefragfile = params.fragments ? file(params.fragments) : null
thecellfile = params.cellcsv   ? file(params.cellcsv)   : null
thebamfile  = params.posbam    ? file(params.posbam)    : null
*/


def just_name(full_file_path) {
  full_file_path.toString().split('/')[-1]           // fixme: use generic path separator
}

process prepare_cr {

  tag "cr-prep $cellbatchsize"

  publishDir "${params.outdir}", pattern: 'cellmetadata',   mode: 'copy'

  when: !params.mermul && !params.muxfile

  input:
  file(f_cells) from ch_cellfile_cr.collect()
  file(posbam) from ch_bamfile_cr.collect()
  val cellbatchsize from params.cellbatchsize
  val winsize from params.winsize

  output:
  file('cellmetadata/cells.tab')    into ch_celltab_cr
  file('cellmetadata/win.tab')      into ch_wintab_cr
  file('cellmetadata/sample.chrlen')    into ch_chrom_length_cr
  set val("crsingle"), file('c_c.*') into ch_demux_cr

  shell:
  filter = params.ncell > 0 ? "head -n ${params.ncell}" : "cat"
  '''
  mkdir -p cellmetadata
# Chrosome length file.
  samtools view -H !{posbam}   \\
    | grep '@SQ'$'\\t''SN:'    \\
    | perl -ne '/\\bSN:(\\S+)/ && ($name=$1); /\\bLN:(\\d+)/ && ($len=$1); print "$name\\t$len\\n";' \\
    | uniq                     \\
    | grep -i 'chr[a-z0-9][a-z0-9]*\\>' \\
    > cellmetadata/sample.chrlen
        # grep -> WARNING DANGERSIGN fixme (see other locations)

# Names + info of selected cells.
  cat !{f_cells} | tr ',' '\\t' | perl -ane 'print if $F[9] == 1' \\
     | (sort -rnk 2 || true) | !{filter} > cellmetadata/chosen_cells.info

# Just the names of selected cells. 
  cut -f 1 cellmetadata/chosen_cells.info | sort > cellmetadata/cell.names

# Tab file for mcx matrix loading
  nl -v0 -nln -w1 < cellmetadata/cell.names > cellmetadata/cells.tab

# Batch lists for demuxing
  split -l !{cellbatchsize} cellmetadata/cell.names c_c.

# Tab file for windows
  ca_make_chromtab.pl !{winsize} cellmetadata/sample.chrlen > cellmetadata/win.tab
  '''
}


process prepare_cr_mux {     // integrate multiple fragment files

  tag "cr-mux-prep $sampleid $sampletag"

  container 'quay.io/cellgeni/cellclusterer'

  when: !params.mermul && !params.posbam

  publishDir "${params.outdir}", pattern: 'cellmetadata',   mode: 'copy'

  input:
  set val(sampletag), val(sampleid), val(root) from ch_mux
        .splitCsv(sep: '\t')
        .view{"manifest: $it"}
  val cellbatchsize from params.cellbatchsize
  val winsize from params.winsize

  output:
  file("cellmetadata/${sampleid}.names")  into ch_cellnames_many_cr
  file('cellmetadata/*-sample.chrlen')    into ch_chrom_length_many_cr
  file('cellmetadata/*-win.tab')          into ch_wintab_many_cr
  file("cellmetadata/${sampleid}.info")   into ch_cellinfo_many_cr
  set val(sampletag), file("*.fragments.gz"), file('c_c.*') into ch_demux_many_cr

  shell:
  filter = params.ncell > 0 ? "head -n ${params.ncell}" : "cat"
  '''
# make fragments, bam, manifest available
  bamfile=!{sampleid}.bam
  fragfile=!{sampleid}.fragments.gz               # output; not used below
  cellfile=!{sampleid}.singlecell.csv

  ln -s !{root}/possorted_bam.bam $bamfile
  ln -s !{root}/fragments.tsv.gz $fragfile
  ln -s !{root}/singlecell.csv $cellfile

  mkdir -p cellmetadata

# Chrosome length file.
  samtools view -H $bamfile  \\
    | grep '@SQ'$'\\t''SN:'    \\
    | perl -ne '/\\bSN:(\\S+)/ && ($name=$1); /\\bLN:(\\d+)/ && ($len=$1); print "$name\\t$len\\n";' \\
    | grep -i 'chr[a-z0-9][a-z0-9]*\\>' \\
    | uniq                    \\
    > cellmetadata/!{sampletag}-sample.chrlen
        # grep -> WARNING DANGERSIGN fixme (see other locations)

# Tab file for windows
  ca_make_chromtab.pl !{winsize} cellmetadata/!{sampletag}-sample.chrlen > cellmetadata/!{sampletag}-win.tab

# Names + info of selected cells.
  perl -F, -ane 's/,/\t/g; print "!{sampletag}-$_" if $F[9] == 1' $cellfile \\
     | (sort -rnk 2 || true) | !{filter} > cellmetadata/!{sampleid}.info

# Just the names of selected cells. 
  cut -f 1 cellmetadata/!{sampleid}.info | sort > cellmetadata/!{sampleid}.names

# Batch lists for demuxing
  split -l !{cellbatchsize} cellmetadata/!{sampleid}.names c_c.

  '''
}


process prepare_mm {        // merge multiplets

  tag "mm-prep $cellbatchsize"

  container 'quay.io/cellgeni/cellclusterer'

  publishDir "${params.outdir}", pattern: 'cellmetadata_mm',   mode: 'copy'

  when: params.mermul && !params.muxfile

  input:
  file(f_cells) from ch_cellfile_mm.collect()
  file(posbam) from ch_bamfile_mm.collect()
  val cellbatchsize from params.cellbatchsize
  val winsize from params.winsize

  output:
  file('cellmetadata/cells.tab')    into ch_celltab_mm
  file('cellmetadata/win.tab')      into ch_wintab_mm
  file('cellmetadata/sample.chrlen')    into ch_chrom_length_mm

  set val("mm"), file('fragmints2.tsv.gz'), file('c_c.*') into ch_demux_mm

  shell:
  '''
  mkdir -p cellmetadata
# Chrosome length file.
  samtools view -H !{posbam}  \\
    | grep '@SQ'$'\\t''SN:'    \\
    | perl -ne '/\\bSN:(\\S+)/ && ($name=$1); /\\bLN:(\\d+)/ && ($len=$1); print "$name\\t$len\\n";' \\
    | uniq | sort -V           \\
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
  cut -f 2 cellmetadata/bc-mm.map | sort -u > cellmetadata/cell.names

# Tab file for mcx matrix loading
  nl -v0 -nln -w1 < cellmetadata/cell.names > cellmetadata/cells.tab

# Batch lists for demuxing
  split -l !{cellbatchsize} cellmetadata/cell.names c_c.

# Tab file for windows
  ca_make_chromtab.pl !{winsize} cellmetadata/sample.chrlen > cellmetadata/win.tab
  '''
}


process join_muxfiles {

// todo fixme this process is activated somehow in single-sample mode,
// even with a draconian `when: false` directive.
// it then sends a tab file and this leads to disaster.
// Need to understand `when` + process triggering (asked Paolo about that once).
// IIRC the issue then was a child or grandchild popping up unexpectedly.
// In this case join_muxfiles is the (grand)child.

  publishDir "${params.outdir}/cellmetadata", pattern: 'singlecell.tsv',   mode: 'link'

  when: false

  input:
  file(fnames) from ch_cellnames_many_cr.toSortedList { just_name(it) }
  file(fninfo) from ch_cellinfo_many_cr.toSortedList { just_name(it) }

  output:
  file('merged.tab') into ch_celltab_manymerged_cr
  file('singlecell.tsv') into ch_seurat_singlecellcsv

  shell:
  '''
  cat !{fnames} > merged.names
  nl -v0 -nln -w1 < merged.names > merged.tab
true
  echo -e "barcode\\ttotal\\tduplicate\\tchimeric\\tunmapped\\tlowmapq\\tmitochondrial\\tpassed_filters\\tcell_id\\tis__cell_barcode\\tTSS_fragments\\tDNase_sensitive_region_fragments\\tenhancer_region_fragments\\tpromoter_region_fragments\\ton_target_fragments\\tblacklist_region_fragments\\tpeak_region_fragments\\tpeak_region_cutsites" > singlecell.tsv
  cat !{fninfo} >> singlecell.tsv
  '''
}


/*
   The merge-multiplet code was inserted in a clone of the preparation process.
   For further excitement, we also cater for merging multiple 10x experiments.
   below are all the merging Y-junctions for plumbing the channels to the same destination.
      prepare_cr      -  pure cr code
      prepare_mm      -  multiplet merging code
      prepare_cr_mux  -  cr code supporting multiple experiments (not yet multipletted)
*/

  ch_celltab_manymerged_cr
    .mix(ch_celltab_cr, ch_celltab_mm)
    .into { ch_celltab; ch_celltab2; ch_celltab3; ch_celltab4; ch_celltab5 }

  ch_wintab_many_cr.toSortedList{ just_name(it) }
    .flatten().first()
    .mix(ch_wintab_cr, ch_wintab_mm)
    .into { ch_wintab; ch_wintab2; ch_wintab3; ch_wintab4 }

  ch_chrom_length_many_cr.toSortedList { just_name(it) }.flatten().first()
    .mix(ch_chrom_length_cr, ch_chrom_length_mm)
    .into { ch_chrom_length; ch_chrom_length2; ch_chrom_length3 }

  ch_fragfile_cr.map { fragf -> [ "crsingle", fragf ] }
    .join(ch_demux_cr)
    .mix(ch_demux_mm, ch_demux_many_cr)           // fixme this block is a bit convoluted
    .transpose()
    .set { ch_demux_batch }


process sample_demux {

  tag "sample-$sampletag batch-$batchtag"

  input:
  set val(sampletag), file(frags), file(cells) from ch_demux_batch.view{"sample_demux: $it"}

  output:
  file('*celldata/[ACGT][ACGT][ACGT][ACGT]/*.bed') into ch_demuxed

  file('*celldata/mtx*.edges') into ch_all_edges
  set val(sampletag), file('*celldata/mtx*.edges') into ch_sample_edges

  shell:
  batchtag = cells.toString() - 'c_c.'
  thesampletag = sampletag
  '''
  dir=celldata
  if [[ -n "!{thesampletag}" ]]; then
    dir=!{thesampletag}-celldata
  fi
  mkdir -p $dir
  zcat !{frags} | samdemux.pl --barcodefile=!{cells} --outdir=$dir --bucket --fragments --ntest=0 --fnedges=mtx.!{thesampletag}-!{batchtag}.edges --tag=!{thesampletag}
  '''
}


                  // fixme noteme; this depends on bed suffix;
                  // does samdemux outputs sam suffix in possorted bam mode.
ch_demuxed
  .flatMap()
  .map { it -> it.toString() - ~/.*celldata\/[ACGT]{4}\// - ~/\.bed/ + '\t' + it.toString()  }
  .into { ch_cellpaths_clusterindex; ch_cellpaths_masterpeakcov }
  

process make_sample_matrix {

  tag "$sampleid"

  container 'quay.io/cellgeni/cellclusterer'
  
  output:
  file("${sampleid}.w2c.mcx") into ch_sample_join_matrix
  file("${sampleid}.stats") into ch_sample_join_window

  input:
  set val(sampleid), file(fnedges) from ch_sample_edges.groupTuple(sort: true).view()
  file('allcelltab') from ch_celltab4.collect()
  file('wintab') from ch_wintab3.collect()

            // fixme: write w2c directly (swap fields)
  shell:
  '''
  export MCLXIOFORMAT=8
	mtx=!{sampleid}.c2w.mcx
  mtxtp=!{sampleid}.w2c.mcx
  cat !{fnedges} | mcxload --stream-split -abc - -strict-tabc allcelltab -restrict-tabr wintab --write-binary -o $mtx
  mcxi /$mtx lm tp /$mtxtp wm
  mcx query -imx $mtxtp -tab wintab | tail -n +2 > !{sampleid}.stats
  '''
}


process join_sample_matrix {

  container 'quay.io/cellgeni/cellclusterer'
 
  input:
  file(w2c) from ch_sample_join_matrix.toSortedList { just_name(it) }
  file(win) from ch_sample_join_window.toSortedList { just_name(it) }
  file(wintab) from ch_wintab4.collect()
  file(celltab) from ch_celltab5.collect()
  val ntfs      from  params.ntfs
  val ntfs_ofs  from  params.ntfs_ofs

  output:
  file('mmtx') into ch_load_mmtx2

/*
Need to create these outputs, as they are currently hardcoded in seurat script (fixme).
    f_binary_mat
    regions.names
    cell.names
*/

  shell:
  '''
  ## ca_winsect rank-transforms the data.
  ca_winsect.pl -1 !{wintab} !{win} > __winsel.stats

	# All approaches: use ranks to mitigate sample-level expression differences.

			# first approach: sum ranks, take top N lowest sumrank windows.
			# This implies taking windows with high expression across samples.
  if !{params.sumrank}; then
    ( set +o pipefail; cut -f 1,2 __winsel.stats | tail -n +!{ntfs_ofs} | head -n !{ntfs} | sort -nk 1 ) > window.tab
  else
		# Second approach try to find windows that are variably expressed.
		# We use IQR approach (on the ranks) skewed towards the higher expressed (== lower rank) windows.
		# This is experimental code. Will be either removed or made clean.
  R --no-save <<EOC
	wintab <- read.table("!{wintab}", as.is=T, colClasses=c("character", "character"))

	N <- !{ntfs}
	rd <- read.table("winsect.stats", header=T, row.names=1, as.is=T)				# rank data
	sumrank <- rd[, "Sumrank"]
	names(sumrank) <- rownames(rd)
	srk_names <- names(sort(sumrank))[1:N]
	rd2 <- rd[,-1]			    # remove sumrank column
	rd2[rd2 > 1*N] <- 1*N		# truncate counts at 1N; to avoid high-rank windows dominating; point of interest/debate
													# Given skewed quantile this probably does not have any impact.
	win_choice <- apply(rd2, 1, function (x) { quantile(x, 0.6) - quantile(x, 0.1) })

	iqr2_names <- names(sort(win_choice, decreasing=TRUE))[1:N]
	sum(srk_names %in% iqr2_names)
	write.table(wintab[ wintab[,"V2"] %in% iqr2_names, ], "window.tab", quote=F, row.names=F, col.names=F, sep="\t")
  write.table(as.data.frame(win_choice), "window.quant", quote=F, row.names=T, col.names=F, sep="\t")
	EOC
  fi

  for m in !{w2c}; do
    # input has format <tag>.w2c.mcx
    subfile=__${m%.w2c.mcx}.sub.mcx
    mcxsubs -imx $m --from-disk -tab window.tab "dom(c, t()), out($subfile,wb)"
    echo "done $m"
  done
  export MCLXIOFORMAT=8

		# below constructs an mcxi command to add up a bunch of matrices.
		# The construction looks horrible, the redeeming feature is that it is really fast
    # and the constructed file is pretty simple.
  (
	echo 'mcxi <<EOC'
	echo __*.sub.mcx | perl -ane '@G = map { "/$_ lm add\\n" } @F; $G[0] =~ s/ add//; print @G;'
  echo tp /c2w.mcx wm pop /w2c.mcx wm
  echo EOC
  ) > makematrix.sh
 	bash -e makematrix.sh
  # ri: re-indexed.
  mcxmap -imx w2c.mcx -make-mapc win.map -o w2c.ri.mcx
  mcxmap -tab window.tab -map win.map -o win.ri.tab
  cut -f 2 win.ri.tab > win.names
  cut -f 2 !{celltab} > cell.names

 >&2 echo "Producing matrixmarket format"
	n_entries=$(mcx query -imx w2c.ri.mcx | tail -n +2 | datamash sum 2)

  mkdir mmtx
  ca_make_mmtx.sh -r win.names -c cell.names -m w2c.ri.mcx \\
      -e $n_entries -t pattern -o mmtx/filtered_window_bc_matrix.mmtx.gz
  ln win.names    filtered_win.txt
  ln win.names    mmtx/regions.names
  ln cell.names   filtered_bc.txt
  ln cell.names   mmtx/cell.names
  '''
}


process make_big_matrix {

  container 'quay.io/cellgeni/cellclusterer'

  publishDir "${params.outdir}/cellmetadata", mode: 'link', pattern: 'cell2win.mcx'

  input:
  file all_edges from ch_all_edges.toSortedList{ just_name(it) }
  file celltab from ch_celltab
  file wintab  from ch_wintab

  output:
  file 'cell2win.mcx' into ch_cell2win
  set file('cell2win.mcx'), file(celltab), file(wintab) into ch_dump_big_matrix

  shell:
  '''
  cat !{all_edges} | mcxload --stream-split -abc - -strict-tabc !{celltab} -restrict-tabr !{wintab} --write-binary -o cell2win.mcx
  sleep 3
  '''
}


process mmtx_big_matrix  {

  tag "raw_window_bc_matrix"

  container 'quay.io/cellgeni/cellclusterer'

  publishDir "${params.outdir}/win_matrix", mode: 'link'

  when: !params.stopatcluster

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


process filter_big_matrix {

  tag "bed-mcx-mmtx"

  container = 'quay.io/cellgeni/cellclusterer'

  publishDir "$params.outdir/win_matrix", mode: 'link', pattern: 'other_publish/filtered*',
    saveAs: { fname -> fname - ~/other_publish\// }

  when: false

  input:
  file('cells.tab') from ch_celltab2
  file('win.tab') from ch_wintab2

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
  ca_top_region.sh            \\
      -c ../cells.tab         \\
      -w ../win.tab           \\
      -m ../cell2win.mcx              \\
      -n !{ntfs}                      \\
      -C ../outputs/filtered_cell.stats          \\
      -W ../outputs/win.stats                    \\
      -R ../outputs/regions.names                \\
      -X ../outputs/filtered_window_bc_matrix.mmtx.gz  \\
      -N ../outputs/cell.names

  ln ../outputs/regions.names    ../other_publish/filtered_win.txt
  ln ../outputs/cell.names      ../other_publish/filtered_bc.txt
  ln ../outputs/filtered_window_bc_matrix.mmtx.gz   ../other_publish
  ln ../outputs/filtered_cell.stats                 ../other_publish
  '''
}


process seurat_clustering {

  tag "seurat2020"

  container = 'quay.io/cellgeni/cellclusterer'

  publishDir "$params.outdir/qc", mode: 'link'

  when: params.usecls == '__seurat__'

  input:
  val nclades   from  params.nclades
  val sampleid  from  params.sampleid
  val npcs      from  params.npcs
  file('singlecell.csv') from ch_cellfile_seurat.mix(ch_seurat_singlecellcsv).collect()
  file('mmtx') from ch_load_mmtx2

  output:
  file('seurat-clades.tsv') into ch_seurat_clades
  file('seurat.pdf')
  file('seurat*.rds')

  shell:
  '''
  R --no-save < !{baseDir}/bin/ca_seurat_clades.R
  #ln -s !{baseDir}/bin/cusanovich2018_lib.r .
  #R --slave --quiet --no-save --args  \\
  #--nclades=!{nclades}                \\
  #--npcs=!{npcs}                      \\
  #--matrix=inputs/filtered_window_bc_matrix.mmtx.gz  \\
  #--regionnames=inputs/regions.names  \\
  #--cellnames=inputs/cell.names      \\
  #--winstats=inputs/win.stats         \\
  #--cellstats=inputs/filtered_cell.stats  \\
  #--sampleid=!{sampleid}              \\
  #< !{baseDir}/bin/cluster_cells_cusanovich2018.R
  '''
}


process cusanovich_clustering {

  tag "cusanovich2018"

  container = 'quay.io/cellgeni/cellclusterer'

  publishDir "$params.outdir/qc", mode: 'link', pattern: 'cus.*'

  when: params.usecls == '__cusanovich__'

  input:
  val nclades   from  params.nclades
  val sampleid  from  params.sampleid
  val npcs      from  params.npcs
  file('inputs') from ch_load_mmtx

  output:
  file('cus.obj.*.clades.tsv') into ch_cusanovich_clades
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
  --cellnames=inputs/cell.names       \\
  --winstats=inputs/win.stats         \\
  --cellstats=inputs/filtered_cell.stats  \\
  --sampleid=!{sampleid}              \\
  < !{baseDir}/bin/cluster_cells_cusanovich2018.R
  '''
}

ch_usercls.mix(ch_cusanovich_clades, ch_seurat_clades).set{ ch_clustering }

process clusters_index {

  tag "${cladefile}"

  when: !params.stopatcluster

  input:
  file metafile from ch_cellpaths_clusterindex.collectFile(name: 'cov.inputs', newLine: true)
  file cladefile from ch_clustering

  output:
  file('clusinfo.cl*') into ch_clusinfo

  shell:
  '''
  splitbyclus.pl !{cladefile} !{metafile} clusinfo.cl
  '''
}

ch_clusinfo.flatMap().map { [ it.baseName - 'clusinfo.cl', it ] }.into { ch_cat_cluster_inputs; ch_per_cluster_inputs }


  // This process derives the cluster tag from the file name. We only need to do this once,
  // in subsequent processes we can pass this tag on.


process clusters_merge_inputs {

  tag "${clustag}"

  publishDir "${params.outdir}/clusdef", mode: 'link'

  input:
  set val(clustag), file(clusmetafile) from ch_cat_cluster_inputs

  output:
  set val(clustag), file("cluster.${clustag}.b??") into ch_clustermacs
  
                      // Fixme noteme b?? captures both bed and bam.  A bit fragile; but bam is
                      // probably not needed anymore, and a new format unlikely. Still, noted.

  shell:
  '''
  cut -f 2 !{clusmetafile} | tr '\\n' '\\0' > clusmetafile0
  sort -m -k 1,1V -k 2,2n -k 3,3n --files0-from=clusmetafile0 | perl -pe 'chomp;$_.="\\t+\\n";' > cluster.!{clustag}.bed
  '''
}


process clusters_macs2 {

  tag "${clustag}"

  container = 'fooliu/macs2'

  publishDir "${params.outdir}/macs2", mode: 'link'

  when: !params.stopatcluster

  input:
  set val(clustag), file(clusregionfile) from ch_clustermacs

  output:
  set file('*.xls'), file('*.bed')
  set val(clustag), file('*.narrowPeak') into ch_combine_clusterpeaks, ch_per_cluster_analysis

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


process peaks_make_masterlist {

  tag "masterlist"

  container = 'quay.io/cellgeni/cellclusterer'

  publishDir "${params.outdir}/peaks", mode: 'link'

  when: params.mergepeaks

  input:
  file np_files from ch_combine_clusterpeaks.map { it[1] }.toSortedList { just_name(it) }
  file sample_idxstats from ch_chrom_length.collect()

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

  container = 'quay.io/cellgeni/cellclusterer'

  cache 'lenient'

  // publishDir "${params.outdir}/mp_counts"
  // A lot of files. This is simply the raw input for the cell/peak matrix.

  input:
  file(masterbed_sps) from ch_masterbed_sps.collect()
  file sample_chrlen from ch_chrom_length2.collect()

  file(celldef_list) from ch_cellpaths_masterpeakcov
    .toSortedList { just_name(it) }
    .flatMap()
    .collate(params.cellbatchsize)
    .map { it.join('\n') + '\n' }

  output:
  file('*.mp.txt') into ch_cellpeak

  shell:
  '''
  while read cellname celldef; do
celldef2=ttt.$cellname
sort -k 1,1V -k 2,2n -k 3,3,n $celldef > celldef2
    bedtools coverage               \\
    -a !{masterbed_sps}             \\
    -b $celldef2 -sorted -header     \\
    -g !{sample_chrlen} | awk -F"\t" '{if($4>0) print $0}' > $cellname.mp.txt
  done < !{celldef_list}
  '''
}


process make_subset_peakmatrix {

  tag "${clustag} subset-peak-cell-matrix"

  container = 'quay.io/cellgeni/cellclusterer'

  publishDir "${params.outdir}/clus_peak_matrix", mode: 'link'

  when: params.perclusterpeaks

  input:
  set val(clustag), file(clusmetafile), file(npeakfile) from ch_per_cluster_inputs.join(ch_per_cluster_analysis)
  file sample_chrlen from ch_chrom_length3.collect()

  output:
    // file('cell2peak.gz')
  file('*.peaks_bc_matrix.mmtx.gz') optional true
  file('*.bc_peaks_matrix.mmtx.gz') optional true
  file('*.peaks.txt') optional true
  file('*.bc.txt') optional true
  file('*.macs2-empty') optional true

  shell:
  '''
  if [[ ! -s !{npeakfile} ]]; then
    > !{clustag}.macs2-empty
    exit 0
  fi
  cut -f 1-3 !{npeakfile} | sort -k1,1 -k2,2n -k 3,3n > clusterpeak.bed

  bedtools sort -faidx !{sample_chrlen} -i clusterpeak.bed > clusterpeak_sps.bed
          # ^ similar to peaks_makemasterlist; we only need the selection of columns and sorting


  while read cellname celldef; do
    bedtools coverage               \\
    -a clusterpeak_sps.bed          \\
    -b $celldef -sorted -header     \\
    -g !{sample_chrlen} | awk -F"\t" '{if($4>0) print $0}' > $cellname.mp.txt

    echo -e "$cellname\t$cellname.mp.txt" >> peak.inputs
          # construct this file within-process.
  done < !{clusmetafile}
  cut -f 1 peak.inputs > cellnames.txt
          # idem construct this within-process.
          # ^ similar to cells_masterlist_coverage 

  ca_peak_matrix.sh -c cellnames.txt -p clusterpeak_sps.bed -i peak.inputs -x "!{clustag}."
          # ^ similar to make_master_peakmatrix
  '''
}



process make_master_peakmatrix {

  tag "peak-cell-matrix"

  container = 'quay.io/cellgeni/cellclusterer'

  publishDir "${params.outdir}/peak_matrix", mode: 'link'

  when: !params.stopatcluster

  input:
  file metafile from ch_cellpeak.flatMap { ls -> ls.collect{ it.toString() } }.collectFile(name: 'peak.inputs', newLine: true)
  file('masterpeak.bed') from ch_masterbed_sps2.collect()
  file('cells.tab') from ch_celltab3.collect()

  output:
    // file('cell2peak.gz')
  file('peaks_bc_matrix.mmtx.gz')
  file('bc_peaks_matrix.mmtx.gz')
  file('peaks.txt')
  file('bc.txt')

  shell:
  '''
  cut -f 2 cells.tab > cellnames.txt
                  # fixme: hardcoded output names.
  ca_peak_matrix.sh -c cellnames.txt -p masterpeak.bed -i peak.inputs
  '''
}



