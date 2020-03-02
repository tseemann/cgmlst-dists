#!/usr/bin/env perl
use strict;
use Data::Dumper;

my $DEBUG=1;

my @gene;
my @call;
my $row=0;
my @id;

while (<ARGV>) {
  chomp;
  my($id,@row) = split m/\t/;
  push @id, $id;
  push @call, \@row;
  print STDERR "\rRow ",++$row;
#  last if $DEBUG and $row > 100;
}
shift @id; # remove 'FILE' header
my @gene = @{ shift @call }; $row--;
my $col = scalar(@gene);
print STDERR "\rLoaded $row rows x $col columns\n";

my @dist;
for my $j (0 .. $#id) {
  for my $i ($j .. $#id) {
    $dist[$j][$i] = distance( $call[$j], $call[$i] );
    print STDERR "\rCalling $j vs $i => $dist[$j][$i]";
  }
}
print STDERR "\rWriting distance matrix...\n";

print tsv('DISTANCES', @id);
for my $d (@dist) {
  print tsv(shift(@id), @$d);
}

print STDERR "Done.\n";

sub distance {
  my($x,$y) = @_;
  my $diff=0;
  for my $k (0 .. $#$x) {
    #print STDERR "[$k] $x->[$k] vs $y->[$k]\n" if $DEBUG;
    $diff++ if $x->[$k] ne $y->[$k];
  }
  return $diff;
}

sub tsv { join("\t", @_)."\n"; }



__DATA__
FILE    SC0831.fasta    SEN0401.fasta   SNSL254_A1382.fasta     SPAB_04503>
2019-31278      1       3       LNF     2       LNF     LNF     LNF     1 >
2013-13488      1       3       LNF     2       LNF     LNF     LNF     1 >
2003-20819      1       9       LNF     3       LNF     LNF     13      7 >
2018-26645      2       208     LNF     3       LNF     LNF     304     26>
2019-24143-1    2       19      LNF     3       LNF     LNF     23      16>
2019-28117      1       3       LNF     2       LNF     LNF     LNF     1 >
2019-24338      2       19      LNF     3       LNF     LNF     23      16>