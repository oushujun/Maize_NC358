#!/usr/bin/perl -w
use strict;
#Get Ngap length in a region
#Shujun Ou (shujun.ou.1@gmail.com)
#09/30/2019

my $target = $ARGV[0];
my $Ngap = $ARGV[1];
my $anno = $ARGV[2];
my $rm13N = 1; #1 will discard 13N gaps

open TGT, "<$target" or die $!;
while (<TGT>){
	chomp;
	my ($chr, $s, $e)=(split);
	my @gap = `awk '{if (\$1=="$chr" && \$2>=$s && \$3<=$e) print \$0}' $Ngap`;
	my $length = 0;
	foreach (@gap){
		my $len = (split)[3];
		next if $len eq 13 and $rm13N eq 1;
		$length += $len;
		}
	print "$_\t$length\t$anno\n";
	}
close TGT;

#target list:
#chr9 11645817 11934876

#Ngap list:
#chr9    11639509        11639521        13
#chr9    11732358        11867021        134664
