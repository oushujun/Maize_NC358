#!/usr/bin/perl -w
use strict;
#Shujun Ou (05/12/2019, shujun.ou.1@gmail.com)

my $usage = "
	Count uniq and multi-mapped reads in window-basis
	Determine if a window is homozygous or heterozygous

	perl bam_read_density.pl sorted.bam > sorted.bam.density.xls\n\n";
my $file = $ARGV[0];
die $usage unless -s $file;

open File, "samtools view $file|" or die $usage;
my %pool;
my $start = 1;
my $win = 50000;
my $step = 50000;
my $diff = 10; #if uniq - multi mapped read count >= $diff, then classify this as heterozygous, otherwise homozygous
my $uniq_count = 0;
my $dup_count = 0;
my $scaf = "NA";
my $type = "NA";
print "#Scaffold\tStart\tEnd\tUniq_count\tMulti_count\tType\n";
while (<File>){
	my ($id, $pos, $MQ, $cigar) = (split)[2,3,4,5];
	$id = "scaf_$id";
	my $NH = $1 if /NH:i:([0-9]+)/;
	my $end = $start + $win;
	$scaf = $id if $scaf eq "NA";

	#output the last window of a scaffold and switch to the new one
	if ($scaf ne $id){
		if (($dup_count - $uniq_count >= $diff) or ($uniq_count == 0 and $dup_count >= 3)){
			$type = "het";
			} else {
			$type = "homo";
			}
		print "$scaf\t$start\t$end\t$uniq_count\t$dup_count\t$type\n" unless $uniq_count == 0 and $dup_count == 0;
		$scaf = $id;
		$start = 1;
		$end = $start + $win;
		$uniq_count = 0;
		$dup_count = 0;
		}

	#only retain 100% matched reads
#	next unless $MQ == 60;
	next unless $cigar eq "100M";

	#exclude overlapping reads
	if (exists $pool{$pos}){
		next;
		} else {
		$pool{$pos} = $pos;
		}

	#calculate window based read count
	if ($end >= $pos){
		$uniq_count ++ if $NH == 1;
		$dup_count ++ if $NH > 1;
		} else {
		if (($dup_count - $uniq_count >= $diff) or ($uniq_count == 0 and $dup_count >= 3)){
			$type = "het";
			} else {
			$type = "homo";
			}
		print "$id\t$start\t$end\t$uniq_count\t$dup_count\t$type\n" unless $uniq_count == 0 and $dup_count == 0;
		$uniq_count = 0;
		$dup_count = 0;
		$start += $step;
		}
	}
close File;

