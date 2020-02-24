#!/usr/bin/perl

## This script implements one way of finding a consensus set of windows from a panel
## of windows sets, each set corresponding to an individual sample.

## The current method is to compute the ranks of windows on a per-sample basis,
## then do a summed-rank approach. It will be replaced
## when/if something better comes along.

## This approach can be contrasted with that where all cell-window edges from all samples
## are streamed into one giant matrix and the top windows are chosen from that matrix
## (based on how many cells a window) has.

## Anecdotal evidence says that the two approaches lead to a 95% in overlap.



use strict;
use warnings;

my $n_windows_required = shift || die "Need top-frequent site number";
my $fgenometab = shift || die "Need global genome tab file";

my @FILENAMES = @ARGV;

open(STATS, ">winsect.stats") || die "Cannot open stats file";

sub read_key_val {
  my $f = shift;
  open(F, "<$f") || die "No $f";
  return { map { chomp; (split())[0,1] } <F> };
}

sub read_win_stats {
  my $f = shift;
  my $set = read_key_val($f);
  my %out_ranked = ();
  my $i = 1;

  ## Make sure that all sets are fully defined on the full tab file,
  ## so that sum-of-rank makes sense.

  for my $k (sort { $set->{$b} <=> $set->{$a} } keys %$set) {
    $out_ranked{$k} = $i;
    $i++ if $i < $n_windows_required;
    # ^ Beyond Do not punish a window beyond a certain rank, for now set at this.
  }

  close(F);
  return \%out_ranked;
}


my %gtab = reverse %{read_key_val($fgenometab)};

  ## Read sample matrix top windows based on #cells with that window
  ## The stored value is the rank based on the same.

my @sets = ();
my $i = 1;
while (my $f = shift) {
  push @sets, [ $i++, read_win_stats($f) ];
}


  ## Print intersection statistics; how many windows do samples share?
  ## For each window compute the sum of ranks.

my %global_sum_rank = ();

for my $x (@sets) {
  my ($ix, $setx) = @$x;
  for my $k (keys %$setx) {
    $global_sum_rank{$k}{total} += $setx->{$k};
  }
}

my @global_win = sort { $global_sum_rank{$a}{total} <=> $global_sum_rank{$b}{total} } keys %global_sum_rank;
my @global_select = @global_win[0..($n_windows_required-1)];
my $nglobal = @global_win;
my $g0 = $global_select[0];
my $g1 = $global_select[1];
print STDERR "$g0 $global_sum_rank{$g0}{total}\n";
print STDERR "$g1 $global_sum_rank{$g1}{total}\n";

for my $x (@sets) {
  my ($ix, $setx) = @$x;
  my $nshare = grep { defined($setx->{$_}) } @global_select;
  printf STDERR "%8d\n", $nshare;
}

print STDERR "$nglobal total\n";

local $" = "\t";
print STATS "Window\tSumrank\t@FILENAMES\n";

for my $w (@global_select) {
  print "$gtab{$w}\t$w\t$global_sum_rank{$w}{total}\n";
  my @stats = map { $_->[1]{$w} } @sets;
  print STATS "$w\t$global_sum_rank{$w}{total}\t@stats\n";
}

close(STATS);

