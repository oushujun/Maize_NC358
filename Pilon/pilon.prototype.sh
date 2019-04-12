#!/bin/bash -login
#SBATCH -N 1
#SBATCH -t 21-0:00:00
#SBATCH --ntasks-per-node 36
#SBATCH --mem=360GB

threads=36	#num of threads

#path to files and programs
PATH=/home/oushujun/las/TE/MaizeNAM/NC358/B73_falcon_only	#path of the project
FILE=B73_falcon_only_consensus_dft.fasta	#fasta
PE1=/home/oushujun/las/TE/MaizeNAM/NC358/illumina/NC358_S0_R1_001.fastq.gz	#fastq.gz
PE2=/home/oushujun/las/TE/MaizeNAM/NC358/illumina/NC358_S0_R2_001.fastq.gz	#fastq.gz

#path to programs
minimap2=/home/oushujun/las/bin/minimap2/minimap2-2.16_x64-linux/minimap2
samtools=/work/LAS/mhufford-lab/oushujun/bin/miniconda2/envs/mapping/bin/samtools
pilon=/work/LAS/mhufford-lab/oushujun/bin/miniconda2/envs/mapping/bin/pilon
call_seq=/home/oushujun/las/git_bin/Maize_NC358/bin/output_by_list.pl
split_scaf=/home/oushujun/las/git_bin/Maize_NC358/GapFix/make_contig_from_scaffold.pl
grep=/usr/bin/grep
perl=/usr/bin/perl


cd $PATH

#check files and programs
if [ ! -s $FILE ]; then
	echo "Input genome file is not existant!"
	exit
elif [ ! -s $PE1 ] || [ ! -s $PE2 ]; then
	echo "At least one of the input fastq files is not existant!"
	exit
elif [ ! -s $minimap2 ] || [ ! -s $samtools ]; then
	echo "At least one of the dependent programs is not existant!"
	exit
fi

#load environment
conda activate mapping

#split scaffolds into contigs
$perl $split_scaf ${FILE} > ${FILE}.contig
FILE=${FILE}.contig

#get the list of all contigs
$grep \> ${FILE} | $perl -nle 's/>//g; print $_' > ${FILE}.list

#index the genome with minimap2
$minimap2 -d ${FILE}.mmi ${FILE}

#align pair-end reads and sort on per job/lane/sample basis
mkdir ILLUMINA
$minimap2 -ax sr -t $threads ${FILE}.mmi $PE1 $PE2 | $samtools view -bS - | $samtools sort -@ $threads - -T ILLUMINA/${FILE}.temp -o ILLUMINA/${FILE}.bam

#merge bam files for the same genome
#$samtools merge ILLUMINA/${FILE}.merge.bam ILLUMINA/*bam
#$samtools index  ILLUMINA/${FILE}.merge.bam


mkdir CONTIG
cd CONTIG
mkdir OUTPUT


#extract per contig bam and fasta, then perform pilon on each contig
for CONTIG in `cat ../${FILE}.list`; do

echo "#!/bin/bash -login
#SBATCH -N 1
#SBATCH -t 21-0:00:00
#SBATCH --ntasks-per-node 36
#SBATCH --mem=360GB

cd $PATH/CONTIG

#load environment
conda activate mapping

#make bam files per contig
$samtools view -b ../ILLUMINA/${FILE}.bam ${CONTIG} > OUTPUT/${CONTIG}.bam

#index bam files
$samtools index OUTPUT/${CONTIG}.bam

#make fasta files per contig
echo ${CONTIG} | $perl $call_seq 1 ${CONTIG} 1 - -FA > OUTPUT/${CONTIG}.fasta

#index fasta files
$samtools faidx OUTPUT/${CONTIG}.fasta

#create one job per contig to split bam/fasta and run pilon. use 4 cpu per contig and request 8 Gb per core. Most contigs are small and finish in a few minutes, but contig1 may take ~30 minutes or more depending on how big it is.
$pilon --genome OUTPUT/${CONTIG}.fasta --bam OUTPUT/${CONTIG}.bam --fix bases --mindepth 10 --minmq 30 --threads $threads --output OUTPUT/${CONTIG}.pilon --changes --verbose

#polishing report
cat OUTPUT/${CONTIG}.pilon.changes  | sed 's/:/\\t/' | sed 's/-/\\t/' | awk '{print \$1\"\\t\"\$2\"\\t\"\$2+100}' > OUTPUT/${CONTIG}.BED
" > worker.${CONTIG}.qsub

#submit per contig jobs
sbatch worker.${CONTIG}.qsub

done

#make features 100 bp in length starting at change in BED that are readily visible in IGV and concatenate all into single file 
#cat OUTPUT/*BED > ${FILE}.all.BED


