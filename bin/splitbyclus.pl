#!/usr/bin/perl

use strict;
use warnings;

my $cladefile   = shift || die "Need clade file (tsv-format: <barcode><tab><clusid><cr>)";
my $bamlistfile = shift || die "Need file with locations of per-sample bam files";
my $prefix = shift || die "Need prefix for files to be created";

open(CLADES, "<$cladefile") || die "no clades";
open(BAMLIST, "<$bamlistfile") || die "no list";

my %lists = ();

my %tag2file = map { chomp; split } <BAMLIST>;

while (<CLADES>) {
  chomp;
  my ($tag, $clusno) = split;
  die "No file for tag $tag" unless defined($tag2file{$tag});
  push @{$lists{$clusno}}, $tag;
}

for my $no (sort { $a <=> $b } keys %lists) {
  my $batchfile = "$prefix$no.txt";
  open(OUT, ">$batchfile") || die "no batch file";
  print OUT map { "$_\t$tag2file{$_}\n" } @{$lists{$no}};
  close(OUT);
}


