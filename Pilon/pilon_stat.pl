#!/usr/bin/perl -w
use strict;

#Insertion: 1bp, >1bp
#Deletion: 1bp, >1bp
#Substitution: 1bp, >1bp

my $usage = "";

my $file = $ARGV[0];
open File, "<$file" or die $usage;

my %changes = ('Ins_1bp_count', '', 'Ins_1bp_seq', '', 'Ins_multi_count', '', 'Ins_multi_seq', '', 'Del_1bp_count', '', 'Del_1bp_seq', '', 'Del_multi_count', '', 'Del_multi_seq', '', 'Sub_1bp_count', '', 'Sub_1bp_bf_seq', '', 'Sub_1bp_af_seq', '', 'Sub_multi_count', '', 'Sub_multi_bf_seq', '', 'Sub_multi_af_seq', '');
while (<File>){
	my ($before, $after) = (split)[2,3];
	my ($len_bf, $len_af) = (length $before, length $after);
	if ($before eq "."){
		if ($len_af == 1){
			$changes{'Ins_1bp_count'}++;
			$changes{'Ins_1bp_seq'} .= $after;
			} else {
			$changes{'Ins_multi_count'}++;
			$changes{'Ins_multi_seq'} .= $after;
			}
		} 
	elsif ($after eq "."){
		if ($len_bf == 1){
			$changes{'Del_1bp_count'}++;
			$changes{'Del_1bp_seq'} .= $before;
			} else {
			$changes{'Del_multi_count'}++;
			$changes{'Del_multi_seq'} .= $before;
			}
		}
	else {
		if ($len_bf == 1 and $len_af == 1){
			$changes{'Sub_1bp_count'}++;
			$changes{'Sub_1bp_bf_seq'} .= $before;
			$changes{'Sub_1bp_af_seq'} .= $after;
			} else {
			$changes{'Sub_multi_count'}++;
			$changes{'Sub_multi_bf_seq'} .= $before;
			$changes{'Sub_multi_af_seq'} .= $after;
			}
		}
	}
close File;

print "Change_type\t1bp\t%A\t%T\t%C\t%G\t>1bp\t%A\t%T\t%C\t%G\n";
my $Ins_1bp_sum = &stat($changes{'Ins_1bp_seq'});
my $Ins_multi_sum = &stat($changes{'Ins_multi_seq'});
my $Del_1bp_sum = &stat($changes{'Del_1bp_seq'});
my $Del_multi_sum = &stat($changes{'Del_multi_seq'});
my $Sub_1bp_bf_sum = &stat($changes{'Sub_1bp_bf_seq'});
my $Sub_1bp_af_sum = &stat($changes{'Sub_1bp_af_seq'});
my $Sub_multi_bf_sum = &stat($changes{'Sub_multi_bf_seq'});
my $Sub_multi_af_sum = &stat($changes{'Sub_multi_af_seq'});

print "Insertion\t$changes{'Ins_1bp_count'}\t$Ins_1bp_sum\t$changes{'Ins_multi_count'}\t$Ins_multi_sum\n";
print "Deletion\t$changes{'Del_1bp_count'}\t$Del_1bp_sum\t$changes{'Del_multi_count'}\t$Del_multi_sum\n";
print "Substitution_before\t$changes{'Sub_1bp_count'}\t$Sub_1bp_bf_sum\t$changes{'Sub_multi_count'}\t$Sub_multi_bf_sum\n";
print "Substitution_after\t$changes{'Sub_1bp_count'}\t$Sub_1bp_af_sum\t$changes{'Sub_multi_count'}\t$Sub_multi_af_sum\n";

#calculate percent of ATCG in a given sequence
sub stat {
	my $seq = $_[0];
	$seq = '' unless defined $seq;
	my $len = length $seq;
	$len = 1 unless $len != 0;
	my $A = $seq =~ tr/Aa//;
	my $T = $seq =~ tr/Tt//;
	my $C = $seq =~ tr/Cc//;
	my $G = $seq =~ tr/Gg//;
	$A = sprintf("%.2f", $A*100/$len);
	$T = sprintf("%.2f", $T*100/$len);
	$C = sprintf("%.2f", $C*100/$len);
	$G = sprintf("%.2f", $G*100/$len);
	my $sum = "$A\t$T\t$C\t$G";
	return $sum;
	}


