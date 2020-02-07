#!/usr/bin/perl

use strict;
use warnings;

my $winsize = shift || die "Need window size";
die "Need integer window size" unless $winsize + 0 eq $winsize;

my $n = 0;

while (<>) {
  chomp;
  my ($name, $len) = split;
  my $i = 0;
  while ($i < $len + $winsize) {      # Add one overflow window as defensive programming.
                                      # Should be empty always and discarded later.
    my $id = $name . '_' . $i;
    print "$n\t$id\n";
    $i += $winsize;
    $n++;
  }
}

