#!/usr/bin/perl -w
use strict;
#Shujun Ou (shujun.ou.1@gmail.com)
##04/11/2019
#Usage:
#	Run this when scripts made by make_batch_qsub.pl are all finished.
#	make script for each genome by perl post_pilon.pl list
#	Run each script and get the polished genome.

#path to files and programs
my $root = "/home/oushujun/las/TE/MaizeNAM/NC358";
my $stitch_contig = "/home/oushujun/las/git_bin/Maize_NC358/GapFix/make_scaffold_from_contig.pl";
my $assembly_stat = "/home/oushujun/las/git_bin/Maize_NC358/bin/assemblathon_stats.pl";
my $pilon_stat = "/work/LAS/mhufford-lab/oushujun/git_bin/Maize_NC358/Pilon/pilon_stat.pl";


while (<>) {
	chomp;
	next if /^$/;
	my $info = $_;
	my ($PATH, $FILE)=($1, $2) if $info=~/(.*)\/(.*)/;
	my $jobtype="stitch.pilon";
	my $threads=36;
	open Qsub, ">worker.$jobtype.$PATH.$FILE.qsub";
	print Qsub "#!/bin/bash -login
#SBATCH -N 1
#SBATCH -t 21-0:00:00
#SBATCH --ntasks-per-node $threads
#SBATCH --mem=160GB

cd $root/$PATH

#concatenicate pilon polished contigs into a multi-fasta file
rm ${FILE}.contig.pilon ${FILE}.contig.pilon.changes ${FILE}.contig.pilon.BED 2>/dev/null
for CONTIG in `cat ${FILE}.contig.list`; do
	if [ -s CONTIG/OUTPUT/\${CONTIG}.pilon.fasta ]; then
		cat CONTIG/OUTPUT/\${CONTIG}.pilon.fasta | perl -nle 's/_pilon//g; print \$_' >> ${FILE}.contig.pilon
		cat CONTIG/OUTPUT/\${CONTIG}.pilon.changes >> ${FILE}.contig.pilon.changes
		cat CONTIG/OUTPUT/\${CONTIG}.BED >> ${FILE}.contig.pilon.BED
		rm CONTIG/OUTPUT/\${CONTIG}.bam CONTIG/OUTPUT/\${CONTIG}.bam.bai CONTIG/OUTPUT/\${CONTIG}.fasta CONTIG/OUTPUT/\${CONTIG}.fasta.fai 2>/dev/null
	else
		echo \"ERROR: The polished contig CONTIG/OUTPUT/\${CONTIG}.pilon.fasta is not exist or has zero size!\"
	fi
done


#stitch contigs back into scaffolds
perl $stitch_contig ${FILE}.contig.pilon > ${FILE}.pilon

#calculate genome stats before and after polishing
perl $assembly_stat ${FILE} > ${FILE}.stat
perl $assembly_stat ${FILE}.pilon > ${FILE}.pilon.stat

#calculate pilon polishing summary
perl $pilon_stat ${FILE}.contig.pilon.changes > ${FILE}.contig.pilon.changes.stat

";
close Qsub;
	}


