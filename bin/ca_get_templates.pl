#!/usr/bin/perl

# This script filters per-cell bam files, following to some extent 10x processing.
# It takes reads that are properly mapped (bit 2) and paired (bit 1) and have
# quality > 30.
# Optionally it filters out duplicates.
# The 10x software groups duplicates and associates them as a count with
# one chosen representative.
# We have not yet implemented this grouping. Currently the options are
# (1) output 'duplicates' as separate regions.
# (2) do not output 'duplicates'; they will not contribute to depth either.

# 10x atac seq fragment processing, quoted from
# https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/algorithms/overview#dups
#
#   "  A barcoded fragment may get sequenced multiple times due to PCR                    "
#   "  amplification. We mark duplicates in order to identify the original fragments      "
#   "  that constitute the library and contribute to its complexity.                      "
#   "  We find duplicate reads by identifying groups of read-pairs, across all barcodes,  "
#   "  where the 5' ends of both R1 and R2 have identical mapping positions on the        "
#   "  reference, correcting for soft-clipping. These groups of read-pairs arise          "
#   "  from the same original molecule. Among these read-pairs, the most common           "
#   "  barcode sequence is identified. One of the read-pairs with that barcode            "
#   "  sequence is labelled the 'original' and the other read-pairs in the group are      "
#   "  marked as duplicates of that fragment in the BAM file. If it passes the            "
#   "  filters described in the next paragraph, this is the only read pair that is        "
#   "  reported as a fragment in the fragment file, and it will be marked with the        "
#   "  most-common barcode sequence.                                                      "
    
#   "  While processing the group of identically aligned read-pairs as described          "
#   "  above, once the original fragment is marked, we determine if the fragment is       "
#   "  mapped with MAPQ > 30 on both reads, is not mitochondrial, and not                 "
#   "  chimerically mapped. If the fragment passes these filters, we create one           "
#   "  entry in the fragments.tsv.gz file marking the start and end of the fragment       "
#   "  after adjusting the 5' ends of the read-pair to account for transposition,         "
#   "  during which the transposase occupies a region of DNA 9 base pairs long (see       "
#   "  figure). With this entry we associate the most common barcode observed for         "
#   "  the group of read-pairs and the number of times this fragment is observed in       "
#   "  the library (size of the group). Note that as a consequence of this approach,      "
#   "  each unique interval on the genome can be associated with only one barcode.        "
#   "  Each entry is tab-separated and the file is position-sorted and then run           "
#   "  through the SAMtools tabix command with default parameters.                        "


use strict;
use warnings;
use Getopt::Long;

my @ARGV_COPY  = @ARGV;
my $n_args = @ARGV;

my $help  =  0;
my $fh      =  \*STDIN;
my $fn      =  "";
my $progname = "ca-get-templates.pl";
my $take_dups = 1;
my $mark_dups = 0;


sub help {
   print <<EOH;
Usage:
   $progname [options]
Options:
--help                  this
--input <filename>      (default reads from STDIN)
--nodups                filter out reads marked as duplicates by 10x
--markdups              mark duplicates (for manual inspection)
EOH
}

if
(! GetOptions
   (  "help"            =>   \$help
   ,  "input=s"         =>   \$fn
   ,  "dups!"           =>   \$take_dups
   ,  "markdups!"       =>   \$mark_dups
   )
)
   {  print STDERR "option processing failed\n";
      exit(1);
   }
if ($help) {
  help();
  exit(0);
}

if ($fn) {
  open(INPUT, "<$fn") || die "cannot read $fn";
  $fh = \*INPUT;
}


%::bits =
( 1    => 'paired'
, 2    => 'proper-mapped'
, 4    => 'unmapped'
, 8    => 'mate-unmapped'
, 16   => 'RRS'
, 32   => 'MRS'
, 64   => 'first'
, 128  => 'second'
, 256  => 'secondary'
, 512  => 'fail_qual'
, 1024 => 'duplicate'
, 2048 => 'supplementary'
) ;

%::stib = reverse %::bits;
%::df = ();


while (<$fh>) {
  chomp;
  my @F = split "\t";
  if (/^@/) {
    # print;
    next;
  }


  my $name = $F[0];
  my $flag = $F[1]+0;
  my $qual = $F[4]+0;

  if (  ( $flag & $::stib{'proper-mapped'}  )
    &&  ( $flag & $::stib{'paired'}         )
    && ( $take_dups || !($flag & $::stib{'duplicate'}) )
    &&  $qual > 30
  ) {
    push @{$::df{$name}}, [ @F ];
  }
}


$" = "\t";

for my $name (sort keys %::df) {
  my @lines = @{$::df{$name}};
  my $n     = @lines;
  if ($n != 2) {
    print STDERR "$name\t$n entries\n";
  }
  else {
    # print "@$_\n" for @lines;
    my $p1 = $lines[0][3];
    my $p2 = $lines[1][3];
    my $ref = $lines[0][2];
    if ($p1 > $p2) {
      ($p1, $p2) = ($p2, $p1);
    }
    my $q1 = $p1 + 3;
    my $q2 = $p2 + 50 - 6;
    my $dist = $q2 - $q1;
    my $dupfield = "";
    $dupfield = "\t$dist *" if $mark_dups && $lines[0][1] & $::stib{'duplicate'};
    print "$ref\t$q1\t$q2$dupfield\n";
  }
}


# Reminder of sam format:
# queryname     (read)
# flag
# refname       (chromo)
# pos
# mapq          (integer)
# cigar
# rnext         reference name of the mate/next read
# pnext         position of the mate/next read
# length        observed template length
# seq           sequence
# qual          quality
#

