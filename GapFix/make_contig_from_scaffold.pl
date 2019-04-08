#!/usr/bin/perl -w
use strict;
#Shujun Ou (shujun.ou.1@gmail.com)
#04/08/2019

my $usage = "\nSplit scaffolds into contigs based on a preset length of N-string.
\n\tperl make_contig_from_scaffold.pl scaffold.fasta > contig.fasta\n\n";

my $genome = $ARGV[0];
my $min_gap_len = 20; #the minimum length of gaps to make a split

open Genome, "<$genome" or die $usage;

$/ = "\n>";
my $id = 0; #id of contig
while (<Genome>){
	s/>//g;
	my ($name, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	while ($seq =~ /(N{$min_gap_len,})/i){
		my $gap = $1;
		my $gap_len = length $gap;
		my $contig = substr ($seq, 0, index ($seq, $gap));
		$seq = substr ($seq, index ($seq, $gap));
		$seq =~ s/^N+//i;
		print ">$name-ctg$id-N$gap_len\n$contig\n";
		$id++;
		}
	print ">$name-ctg$id-N0\n$seq\n";
	}
close Genome;

