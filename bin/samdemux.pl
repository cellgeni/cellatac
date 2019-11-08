#!/usr/bin/perl

# Code adapted and fleshed out by svd from a magnificent one-page one-liner by
# David Jackson.  Output buffering was added, and output for now changed from
# bam to sam to prevent many processes spawning. Output files are optionally
# distributed across buckets constructed from barcode tags to prevent large
# directories.  Not sure the latter is really useful, but at least 'ls' may
# not ruin your terminal on a rainy day.


use strict;
use warnings;
use autodie;
use Getopt::Long;

# die "Include filtering from 10xINE.pl (based on quality and region) before using this code";

my %files;           # hashes to file handles, handle is barcode.
my %cache;           # buffer reads; flush once at certain size.
my %count;           # count reads per barcode currently in cache.

my $SIZE = 1000;     # cache size
my $N_READ = 0;
my $N_WRITTEN = 0;

my @ARGV_COPY  = @ARGV;
my $n_args = @ARGV;

my $help  =  0;
my $progname = 'samdemux.pl';
my $barcodefile = "";
my $od = "demux-out";
my $ntest = 0;
my $bucket = 0;


sub help {
   print <<EOH;
Usage:
   $progname [options]
Options:
--help                  this
--barcodefile=<FNAME>   file with barcodes
--dir=<DIR>             directory
--ntest=<NUM>           test on NUM reads (e.g. 1000000)
--bucket                distribute the files across directory buckets
EOH
}


if
(! GetOptions
   (  "help"            =>   \$help
   ,  "ntest=i"         =>   \$ntest
   ,  "barcodefile=s"   =>   \$barcodefile
   ,  "outdir=s"        =>   \$od
   ,  "bucket"          =>   \$bucket
   )
)
   {  print STDERR "option processing failed\n";
      help();
      exit(1);
   }

if (!$n_args || $help) {
   help();
   exit(0);
}


sub tag_barcode {
  my $bc = shift;
  my $tag = join "", (split "", $bc)[2,4,6,8];
  if (length($tag) != 4 || $tag !~ /^[ACGT]{4}$/) {
    die "Weird tag [$tag] for barcode [$bc]";
  }
  return $tag;
}

sub flush_lines {
  my $bc = shift;
  my $N  = shift;
  my $ar = $cache{$bc};
  my $fh = $files{$bc};
  for my $e (@$ar) {
    print {$fh} $e;
  }
  my $n = @$ar;
  warn "Discrepancy $n $N\n" unless $n == $N;
  $N_WRITTEN += @$ar;
  $cache{$bc} = [];
  $count{$bc} = 0;
}

-e "$od" || mkdir("$od") || die "Could not make directory $od";

open(BARCODES, '<', $barcodefile) || die "Cannot open barcode file $barcodefile";


while (<BARCODES>) {
  chomp;
  my $bc = $_;
  # open(my $f2,"| samtools view -u - |bamstreamingmarkduplicates level=0 tag=CB | samtools view -b - > $od/$bc.bam");

  my $tag = tag_barcode($bc);
  my $dir = $bucket ? "$od/$tag" : $od;
  -e $dir || mkdir("$dir") || die "Could not make directory $dir";
  my $bucket = "$dir/$bc.sam";

  die "File already opened for barcode $bc" if defined($files{$bc});

  open(my $bcdest,"> $bucket") || die "Could not open barcode bucket $bucket";
  $files{$bc} = $bcdest;
  $cache{$bc} = [];
  $count{$bc} = 0;
                     # ^ open pipes for writing, each associated with a barcode.
}
print STDERR "Done reading barcodes\n";
close BARCODES;


while (<>) {
                         # write headers to each barcode file.
  if(/^@/){
    foreach my $fh (values %files) {
      print {$fh} $_;
    }
  }
                         # Get barcode, write this line to the barcode bucket.
                         # CB:Z cell identifier (Z indicates string type).
                         # Cell ranger documentation:
                         # CR Cell barcode
                         # CY Cell barcode read quality
                         # CB Cell barcode that is error-corrected and
                         #    confirmed against a list of known-good barcode sequences 

  elsif (m{\tCB:Z:(\S+)\b}) {

    my $bc = $1;
    my $cache = $cache{$bc} || next;        # ignore filtered barcodes.
    $N_READ++;

    push @$cache, $_;

    if (++$count{$bc} >= $SIZE) {
      flush_lines($bc, $count{$bc});
    }

    if ($N_READ % 10000 == 0) {             # 100 dots per line, each dot 10,000 reads.
      print STDERR '.';
    }
    if ($N_READ % 1000000 == 0) {           # 1M reads per line.
      printf STDERR " %3d\n", $N_READ / 1000000;
    }

    last if $ntest && $N_READ >= $ntest;
  }
}

                       # Wrap-up. Flush buffers and close barcode files.
for my $bc (sort keys %cache) {
  flush_lines($bc, $count{$bc});
  warn "Issue closing sam file for $bc [$!]" if ! close($files{$bc});
}
warn "closed SAMs\n";
warn "Read $N_READ wrote $N_WRITTEN\n";
warn "Regrettably these numbers are not the same\n" if $N_READ != $N_WRITTEN;

exit 0;


__DATA__
# The original one-liner, with some space added.
perl -nle 'use strict; use autodie; our %h;
    BEGIN{open(my$fh,q(<),shift@ARGV); my$od=shift@ARGV; $od//=q();
      while(<$fh>){chomp; open(my$f2,"| samtools view -u - |bamstreamingmarkduplicates level=0 tag=CB | samtools view -b - > $od/$_.bam");$h{$_}=$f2; }
      close $fh
    }

    if(/^@/){foreach my$fh (values %h){print {$fh} $_ }}
    elsif(m{\tCB:Z:(\S+)\b}){ my$fh=$h{$1}||next; print {$fh} $_;}
    
    END{close $_ foreach values %h; warn "closed BAMs\n"}' ${BARCODES} out/${SAMPLE}
