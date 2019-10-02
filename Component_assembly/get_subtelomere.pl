#!/usr/bin/perl -w
use strict;
#This code was used to identify subtelomere matches in NC358 assemblies using the NC358 21k_75x subtelomere coordinates as guidence.
#Hits will be retained when they fall into the same range (with $fuss) with other basic filters
#Usage: perl get_subtelomere.pl blast.out subtelomere.range.bed > subtelomere.bed
#Shujun Ou (shujun.ou.1@gmail.com) 10/02/2019

my $blast = $ARGV[0];
my $range = $ARGV[1];

my $min_iden = 80;
my $min_len = 80;
my $max_evalue = 1e-5;
my $fuss = 500000; #allow range extend this many bp to account for coordinate variations of diff assemblies

open BLAST, "<$blast" or die $!;
open RANGE, "<$range" or die $!;

my %range;
while (<RANGE>){
	my ($chr, $s, $e) = (split);
	($s, $e) = ($e, $s) if $s > $e;
	push @{$range{$chr}}, [$s, $e];
	}
close RANGE;

while (<BLAST>){
	my ($chr, $iden, $len, $s, $e, $evalue) = (split)[1,2,3,8,9,10];
	next if $iden < $min_iden;
	next if $len < $min_len;
	next if $evalue > $max_evalue;
	($s, $e) = ($e, $s) if $s > $e;
	my $keep = 0;
	if (exists $range{$chr}){
		my @range = @{$range{$chr}};
		foreach (@range){
			my ($start, $end) = @{$_}[0,1];
			($start, $end) = ($start - $fuss, $end + $fuss);
			next if $s < $start or $e > $end;
			$keep = 1;
			}
		}
#	print $_ if $keep == 1;
	print "$chr\t$s\t$e\n" if $keep == 1;
	}


