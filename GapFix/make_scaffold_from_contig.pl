#!/usr/bin/perl -w
use strict;
#Shujun Ou (shujun.ou.1@gmail.com)
#04/08/2019

my $usage = "\nCombine contigs into scaffolds. Contigs must follow the naming criteria in make_contig_from_scaffold.pl.
\n\tperl make_scaffold_from_contig.pl contig.fasta > scaffold.fasta\n\n";

my $genome = $ARGV[0];

open Genome, "<$genome" or die $usage;

$/ = "\n>";
my %genome;
while (<Genome>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	my ($name, $contig_id, $gap_len) = ($1, $2, $3) if $id =~ /^(.*)\-ctg([0-9]+)\-N([0-9]+)$/;
	$seq = $seq."N"x$gap_len;
	$genome{$name}{$contig_id} = $seq;
	}
close Genome;
$/ = "\n";

foreach my $name (keys %genome){
	print ">$name\n";
	foreach my $contig_id (sort {$a<=>$b} keys %{$genome{$name}}){
		print "$genome{$name}{$contig_id}";
		}
	print "\n";
	}

