#!/usr/bin/perl -w
use strict;
use File::Basename;
#Shujun Ou (shujun.ou.1@gmail.com) 04/04/2019
#required blastn, call_seq_by_list.pl

my $usage = "\nPatch gaps an old scaffold with sequences in a new scaffold.
\n\tperl patch_gap.pl old_scaffold.fasta new_scaffold.fasta gap.list
\n\t\tgap.list is a list of gap positions in old_scaffold.fasta with format: scaffoid_id gap_start gap_end
\t\tThe patched file is named old_scaffold.fasta.gapfixed\n\n";

#intputs
my $old_genome = $ARGV[0]; #the error-corrected genome with gap
my $new_genome = $ARGV[1]; #the error-prone genome without gap
my $gap_list = $ARGV[2]; #list of gap positions in $old_genome with format: scaffoid_id gap_start gap_end

#dependencies
my $script_path = dirname(__FILE__);
my $callseq = "$script_path/../bin/call_seq_by_list.pl";
my $blastn = "blastn";

#read gap info
open GAPlist, "<$gap_list" or die $usage;
my %gaplist;
while (<GAPlist>){
	my ($id, $start, $end) = (split);
	my $gap_len = $end - $start + 1;
	$gaplist{$id} = [$start, $end, $gap_len];
	}
close GAPlist;

#call flanking sequence of each gap
`perl -nle 'my (\$id, \$start, \$end)=(split); my (\$left, \$right)=(\$start-1000, \$end+1000); (\$start, \$end)=(\$start+100, \$end-100); my \$list="\${id}_left \$id:\$left..\$start\n\${id}_right \$id:\$end..\$right"; print \$list' $gap_list | perl $callseq - -C $old_genome > $gap_list.flank.fa`;

#find new coordinates of gaps in $new_genome
open Blast, "$blastn -query $gap_list.flank.fa -subject $new_genome -outfmt=6 |" or die $!;
my %patlist;
while (<Blast>){
	my ($query, $subject, $iden, $len, $new_s, $new_e)=(split)[0,1,2,3,8,9];
	my ($id, $old_s, $old_e, $part)=($1,$2,$3,$4) if $query=~/^(.*):([0-9]+)\.\.([0-9]+)\|.*(left|right)$/;
	next unless $id eq $subject; #require align to the same scaffold
	next unless $iden > 95 and $len > 850; #alignment at least 850 bp with at least 95% identity
	next if abs($new_s-$old_s) > 10000000; #old and new coordinates do not differ more than 10 Mb
	$patlist{$id}[0] = $new_e if $part eq "left";
	$patlist{$id}[1] = $new_s if $part eq "right";
	if (defined $patlist{$id}[0] and defined $patlist{$id}[1]){
		my $new_len = $patlist{$id}[1] - $patlist{$id}[0] + 1;
		my $diff = abs($new_len - $gaplist{$id}[2]);
		print "Warning: patch sequence is $diff bp different from the gap!\n" if $diff > 100000;
		}
	}

#patch gaps in $old_genome with sequence from $new_genome
$/ = "\n>";
open Old, "<$old_genome" or die $!;
open New, ">$old_genome.gapfixed" or die $!;
while (<Old>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	if (exists $gaplist{$id}){
		print "Fixing the gap in $id\n\tGap coordinate: $gaplist{$id}[0]..$gaplist{$id}[1]\n\tPatch source: $patlist{$id}[0]..$patlist{$id}[1]\n\n";
		my $seq_gap = `echo '$id $id:$patlist{$id}[0]..$patlist{$id}[1]' | perl $callseq - -C $new_genome`;
		$seq_gap = (split /\n/, $seq_gap, 2)[1];
		$seq_gap =~ s/\s+//g;
		my ($start, $end) = ($gaplist{$id}[0], $gaplist{$id}[1]);
		my $seq_left = substr $seq, 0, $start;
		my $seq_right = substr $seq, $end;
		my $seq_fix = "$seq_left"."$seq_gap"."$seq_right";
		print New ">$id\n$seq_fix\n";
		} else {
		print New ">$id\n$seq\n";
		}
	}
close Old;
close New;

