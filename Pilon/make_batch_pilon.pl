#!/usr/bin/perl -w
use strict;
#Shujun Ou (shujun.ou.1@gmail.com)
#04/11/2019
#Usage:
#	get a list of genomes in format: /path/genome.fa
#	then make script for each genome by perl make_batch_qsub.pl list
#	Run each script; rerun scripts if interrupted, will pick up from the interruption.


#path to files and programs
my $PE1 = "/home/oushujun/las/TE/MaizeNAM/NC358/illumina/NC358_S0_R1_001.fastq.gz";      #fastq.gz
my $PE2 = "/home/oushujun/las/TE/MaizeNAM/NC358/illumina/NC358_S0_R2_001.fastq.gz";      #fastq.gz
my $root = "/home/oushujun/las/TE/MaizeNAM/NC358";

#path to programs
my $minimap2 = "/home/oushujun/las/bin/minimap2/minimap2-2.16_x64-linux/minimap2";
my $samtools = "/work/LAS/mhufford-lab/oushujun/bin/miniconda2/envs/mapping/bin/samtools";
my $pilon = "/work/LAS/mhufford-lab/oushujun/bin/miniconda2/envs/mapping/share/pilon-1.23-0/pilon-1.23.jar";
my $call_seq = "/home/oushujun/las/git_bin/Maize_NC358/bin/output_by_list.pl";
my $split_scaf = "/home/oushujun/las/git_bin/Maize_NC358/GapFix/make_contig_from_scaffold.pl";
my $run_jobs_parallel = "/home/oushujun/las/git_bin/Maize_NC358/bin/run_jobs_parallel.pl";



while (<>) {
	chomp;
	next if /^$/;
	my $info = $_;
	my ($PATH, $FILE)=($1, $2) if $info=~/(.*)\/(.*)/;
	my $jobtype="pilon";
	my $threads=36;
#print "$percent\n";next;
	open Qsub, ">worker.$jobtype.$PATH.$FILE.qsub";
	print Qsub "#!/bin/bash -login
#SBATCH -N 1
#SBATCH -t 21-0:00:00
#SBATCH --ntasks-per-node $threads
#SBATCH --mem=360GB

conda activate mapping

cd $root/$PATH

FILE=${FILE}.contig

#split scaffolds into contigs
if [ ! -s \${FILE}.list ]; then
perl $split_scaf ${FILE} > ${FILE}.contig

#get the list of all contigs
grep \\> \${FILE} | perl -nle 's/>//g; print \$_' > \${FILE}.list
fi

#index the genome with minimap2
if [ ! -s \${FILE}.mmi ]; then
$minimap2 -d \${FILE}.mmi \${FILE}
fi

#align pair-end reads and sort on per job/lane/sample basis
mkdir ILLUMINA
if [ ! -s ILLUMINA/\${FILE}.bam ]; then
$minimap2 -ax sr -t $threads \${FILE}.mmi $PE1 $PE2 | $samtools view -bS - | $samtools sort -@ $threads - -T ILLUMINA/\${FILE}.temp -o ILLUMINA/\${FILE}.bam
fi

#index the bam file
if [ ! -s ILLUMINA/\${FILE}.bam.bai ]; then
$samtools index ILLUMINA/\${FILE}.bam
fi

mkdir CONTIG
cd CONTIG
mkdir OUTPUT

#extract per contig bam and fasta, then perform pilon on each contig
rm worker.list 2>/dev/null
for CONTIG in `cat ../\${FILE}.list`; do

echo \"#!/bin/bash -login
#SBATCH -N 1
#SBATCH -t 21-0:00:00
#SBATCH --ntasks-per-node 4
#SBATCH --mem=32GB

cd $root/$PATH/CONTIG

#load environment
#conda activate mapping

if [ ! -s OUTPUT/\${CONTIG}.pilon.fasta ]; then
##make bam files per contig
$samtools view -b ../ILLUMINA/\${FILE}.bam \${CONTIG} > OUTPUT/\${CONTIG}.bam

##index bam files
$samtools index OUTPUT/\${CONTIG}.bam

##make fasta files per contig
echo \${CONTIG} | perl $call_seq 1 ../\${FILE} 1 - -FA > OUTPUT/\${CONTIG}.fasta

##index fasta files
$samtools faidx OUTPUT/\${CONTIG}.fasta

##create one job per contig to split bam/fasta and run pilon. use 4 cpu per contig and request 8 Gb per core. Most contigs are small and finish in a few minutes, but contig1 may take ~30 minutes or more depending on how big it is.
java -Xms16g -Xmx32g -jar $pilon --genome OUTPUT/\${CONTIG}.fasta --bam OUTPUT/\${CONTIG}.bam --fix bases --mindepth 10 --minmq 30 --threads $threads --output OUTPUT/\${CONTIG}.pilon --changes --verbose

##polishing report
#cat OUTPUT/\${CONTIG}.pilon.changes  | sed 's/:/\\\\t/' | sed 's/-/\\\\t/' | awk '{print \\\$1\\\"\\\\t\\\"\\\$2\\\"\\\\t\\\"\\\$2+100}' > OUTPUT/\${CONTIG}.BED
perl -nle 'my \\\$info = (split)[0]; my (\\\$id, \\\$start) = (\\\$1, \\\$2) if \\\$info =~ /^(.*):([0-9]+)/; my \\\$end = \\\$start+100; print \\\"\\\$id\\\\t\\\$start\\\\t\\\$end\\\"' OUTPUT/\${CONTIG}.pilon.changes > OUTPUT/\${CONTIG}.BED
fi

\" > worker.\${CONTIG}.qsub

echo \"worker.\${CONTIG}.qsub\" >> worker.list

done

#run per-contig jobs in batch
perl $run_jobs_parallel -list worker.list -threads 9

";
close Qsub;
	}


