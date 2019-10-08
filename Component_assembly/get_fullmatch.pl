#!/usr/bin/perl -w
use strict;

my $min_cov = 0.8;
my $min_iden = 0.95;
my $max_dist = 500000; #discard hit if coordinate different from more than this length to the query coordinate
while (<>){
	next if /^\s+?\n/;
	my ($id, $chr, $iden, $len, $start) = (split)[0,1,2,3,8];
	next unless defined $iden;
	next if $iden <= $min_iden;
	my ($id_chr, $id_s, $id_e, $loc) = ($1, $2, $3, $4) if $id =~ /^(.*):([0-9]+)\.\.([0-9]+)\|(.*)$/;
	my $id_len = $id_e - $id_s + 1;
	next if $id_chr ne $chr;
	next if abs($start - $id_s) > $max_dist;
	print "$len/$id_len\t$_" if $len/$id_len >= $min_cov;
	}
