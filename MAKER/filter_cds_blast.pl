#!/usr/bin/perl -w
use strict;
#identify complete query alignment 
#Usage: perl filter_cds_blast.pl query.fa blast.out > query.fa.complete.list
#Shujun Ou (10/09/2019 shujun.ou.1@gmail.com)

my $query = $ARGV[0];
my $blast = $ARGV[1];

my $min_iden = 98;
my $min_cov = 0.8;
my $combine = "~/las/git_bin/EDTA/util/combine_overlap.pl";
my $count = "~/las/git_bin/EDTA/util/count_mask.pl";

open Query, "<$query" or die $!;
open Blast, "<$blast" or die $!;
open Filtered, ">$blast.filtered" or die $!;

$/ = "\n>";
my %query;
while (<Query>){
	s/>//;
	my ($id, $seq) = (split /\n/, $_, 2);
	$id =~ s/\s+.*//g;
	$seq =~ s/\s+//g;
	$query{$id} = length($seq);
	}
$/ = "\n";
close Query;

my %blast;
while (<Blast>){
	chmod;
	my ($id, $chr, $iden, $q_s, $q_e, $s_s, $s_e) = (split)[0,1,2,6,7,8,9];
	next unless defined $iden and $iden =~ /^[0-9.]+$/;
	next if $iden < $min_iden;
	my $id_chr = $1 if $id =~ /\-(chr[0-9]+)\-/;
	next unless defined $id_chr and $id_chr eq $chr;
	print Filtered "$_\n";
	$blast{$id} .= "$id\t$q_s\t$q_e\n";
	}
close Blast;
close Filtered;

foreach my $id (keys %blast){
	my $info = $blast{$id};
	my $exec = "perl $combine <(echo -e \"$info\") temp.blast.tttxx";
	qx(bash -c '$exec' 2> /dev/null) if defined $info;
	my $len = `perl $count temp.blast.tttxx` if -s "temp.blast.tttxx";
	`rm temp.blast.tttxx 2> /dev/null`;
	next unless $len > 0;
	my $full = $query{$id};
	next unless $full > 0;
	my $cov = $len/$full;
	print "$id\t$cov\n" if $cov >= $min_cov;
	}

