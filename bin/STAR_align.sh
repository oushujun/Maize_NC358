#!/bin/bash

#Get STAR binary
export PATH=$PATH:/sonas-hs/ware/hpc/data/program_binaries/STAR-2.5.2b/bin/Linux_x86_64_static/

#wrapper script for STAR star_index
#reference="NC358_21k_75x_GG-mapped.fasta"
#annotation=""
seqmode="paired"
seqFolder="/mnt/grid/ware/hpc_norepl/data/data/kapeel/NAM/RNA-seq/reads/nc358/unmapped"
#compression_type="gzipped"
mkdir star_index STAR_output
#mkdir STAR_output
#outFilterType="${outFilterType}"
outFilterType="BySJout"
#outFilterMultimapNmax="${outFilterMultimapNmax}"
outFilterMultimapNmax=20
#alignSJoverhangMin="${alignSJoverhangMin}"
alignSJoverhangMin=8
#alignSJDBoverhangMin="${alignSJDBoverhangMin}"
alignSJDBoverhangMin=1
#outFilterMismatchNmax="${outFilterMismatchNmax}"
outFilterMismatchNmax=999
#alignIntronMin="${alignIntronMin}"
alignIntronMin=20
#alignIntronMax="${alignIntronMax}"
alignIntronMax=1000000
#alignMatesGapMax="${alignMatesGapMax}"
alignMatesGapMax=1000000
runthis=''

#if [ -n "${annotation}" ]; then runthis="$runthis --sjdbGTFfile ${annotation}"; fi
if [ -n "${outFilterType}" ]; then runthis="$runthis --outFilterType ${outFilterType}"; fi
if [ -n "${outFilterMultimapNmax}" ]; then runthis="$runthis --outFilterMultimapNmax ${outFilterMultimapNmax}"; fi
if [ -n "${alignSJoverhangMin}" ]; then runthis="$runthis --alignSJoverhangMin ${alignSJoverhangMin}"; fi
if [ -n "${alignSJDBoverhangMin}" ]; then runthis="$runthis --alignSJDBoverhangMin ${alignSJDBoverhangMin}"; fi
if [ -n "${outFilterMismatchNmax}" ]; then runthis="$runthis --outFilterMismatchNmax ${outFilterMismatchNmax}"; fi
if [ -n "${alignIntronMin}" ]; then runthis="$runthis --alignIntronMin ${alignIntronMin}"; fi
if [ -n "${alignIntronMax}" ]; then runthis="$runthis --alignIntronMax ${alignIntronMax}"; fi
if [ -n "${alignMatesGapMax}" ]; then runthis="$runthis --alignMatesGapMax ${alignMatesGapMax}"; fi

#STAR --runThreadN 16 --runMode genomeGenerate  --genomeDir star_index --genomeFastaFiles "${reference}"

Input=$(basename ${seqFolder});

if [ $seqmode = paired ];
then
	echo "Paired end mode"

	for x in $seqFolder/*
	do
		if [[ "$x" =~ .*1\.fq\.gz$ ]];
		then
		z=$(basename $x 1.fq.gz)
			echo "file of this interation $z"1.fq.gz" $z"2.fq.gz""
			STAR $runthis --runThreadN 16 --genomeDir star_index --readFilesIn "$Input"/"$z"'1.fq.gz' "$Input"/"$z"'2.fq.gz' --readFilesCommand gunzip -c --outFileNamePrefix "$z"  --outSAMtype BAM SortedByCoordinate --outSAMattrIHstart 0 --outSAMstrandField intronMotif
			mkdir STAR_output/"$z"
                        mv "$z"* STAR_output/"$z"
		fi
		if [[ "$x" =~ .*1\.bz2$ ]];
                then
		z=$(basename $x 1.bz2)
			echo "file of this interation $z"1.bz2" $z"2.bz2""
		        STAR $runthis --runThreadN 4 --genomeDir star_index --readFilesIn "$Input"/"$z"'1.bz2' "$Input"/"$z"'2.bz2' --readFilesCommand bunzip2 -c --outFileNamePrefix "$z" --outSAMtype BAM SortedByCoordinate   --outSAMattrIHstart 0 --outSAMstrandField intronMotif
                        mkdir STAR_output/"$z"
                        mv "$z"* STAR_output/"$z"
                fi
		if [[ "$x" =~ .*1\.fastq$ ]];
		then
		z=$(basename $x 1.fastq)
		        echo "file of this interation $z"1.fastq" $z"2.fastq""
			STAR $runthis --runThreadN 4 --genomeDir star_index --readFilesIn "$Input"/"$z"'1.fastq' "$Input"/"$z"'2.fastq' --outFileNamePrefix "$z" --outSAMtype BAM SortedByCoordinate   --outSAMattrIHstart 0 --outSAMstrandField intronMotif
			mkdir STAR_output/"$z"
                        mv "$z"* STAR_output/"$z"
		fi
		if [[ "$x" =~ .*1\.fq$ ]];
                then
		z=$(basename $x 1.fq)
		       echo "file of this interation $z"1.fq" $z"2.fq""
		       STAR $runthis --runThreadN 4 --genomeDir star_index --readFilesIn "$Input"/"$z"'1.fq' "$Input"/"$z"'2.fq' --outFileNamePrefix "$z" --outSAMtype BAM SortedByCoordinate   --outSAMattrIHstart 0 --outSAMstrandField intronMotif
		       mkdir STAR_output/"$z"
                       mv "$z"* STAR_output/"$z"
                fi
		if [[ "$x" =~ .*1\.fastq.gz$ ]];
                then
		z=$(basename $x 1.fastq.gz)
		       echo "file of this interation $z"1.fastq.gz" $z"2.fastq.gz""
		       STAR $runthis --runThreadN 4 --genomeDir star_index --readFilesIn "$seqFolder"/"$z"'1.fastq.gz' "$seqFolder"/"$z"'2.fastq.gz' --readFilesCommand gunzip -c --outFileNamePrefix "$z" --outSAMtype BAM SortedByCoordinate   --outSAMattrIHstart 0 --outSAMstrandField intronMotif
                       mkdir STAR_output/"$z"
                       mv "$z"* STAR_output/"$z"
                fi
		if [[ "$x" =~ .*1\.fastq.bz2$ ]];
                then
		z=$(basename $x 1.fastq.bz2)
			echo "file of this interation $z"1.fastq.bz2" $z"2.fastq.bz2""
		        STAR $runthis --runThreadN 4 --genomeDir star_index --readFilesIn "$Input"/"$z"'1.fastq.bz2' "$Input"/"$z"'2.fastq.bz2' --readFilesCommand bunzip2 -c --outFileNamePrefix "$z" --outSAMtype BAM SortedByCoordinate   --outSAMattrIHstart 0 --outSAMstrandField intronMotif
                       mkdir STAR_output/"$z"
                       mv "$z"* STAR_output/"$z" 
                fi
		if [[ "$x" =~ .*1\.fq.gz$ ]];
                then
		z=$(basename $x 1.fq.gz)
		        echo "file of this interation $z"1.fq.gz" $z"2.fq.gz""
			STAR $runthis --runThreadN 16 --genomeDir star_index --readFilesIn "$seqFolder"/"$z"'1.fq.gz' "$seqFolder"/"$z"'2.fq.gz' --readFilesCommand gunzip -c --outFileNamePrefix "$z"  --outSAMtype BAM SortedByCoordinate --outSAMattrIHstart 0 --outSAMstrandField intronMotif
                       mkdir STAR_output/"$z"
                       mv "$z"* STAR_output/"$z"
                fi
		if [[ "$x" =~ .*1\.fq.bz2$ ]];
                then
		z=$(basename $x 1.fq.bz2)
		        echo "file of this interation $z"1.fq.bz2" $z"2.fq.bz2""
			STAR $runthis --runThreadN 4 --genomeDir star_index --readFilesIn "$Input"/"$z"'1.fq.bz2' "$Input"/"$z"'2.fq.bz2' --readFilesCommand bunzip2 -c --outFileNamePrefix "$z" --outSAMtype BAM SortedByCoordinate   --outSAMattrIHstart 0 --outSAMstrandField intronMotif
                mkdir STAR_output/"$z"
                mv "$z"* STAR_output/"$z"
		fi
	done
	mkdir STAR_bam
        rsync -av --progress STAR_output/*/*bam STAR_bam
fi

if [ $seqmode = single ];
then
	echo "single end mode"
	
	for x in $Input/*
	do
	if [[ "$x" =~ .*\.bz2$ ]];
        then 
	z=$(basename $x .bz2)
		 echo "file of this interation $x"
		 STAR $runthis --runThreadN 4 --genomeDir star_index --readFilesIn "$x" --readFilesCommand bunzip2 -c --outFileNamePrefix "$z"  --outSAMtype BAM SortedByCoordinate   --outSAMattrIHstart 0 --outSAMstrandField intronMotif
	mkdir STAR_output/"$z"
        mv "$z"* STAR_output/"$z"
	fi
	
	if [[ "$x" =~ .*\.gz$ ]];
        then    
        z=$(basename $x .gz)
                 echo "file of this interation $x"
                 STAR $runthis --runThreadN 4 --genomeDir star_index --readFilesIn "$x" --readFilesCommand gunzip -c --outFileNamePrefix "$z"  --outSAMtype BAM SortedByCoordinate   --outSAMattrIHstart 0 --outSAMstrandField intronMotif
        mkdir STAR_output/"$z"
        mv "$z"* STAR_output/"$z"
	fi
	done
	mkdir STAR_bam
	rsync -av --progress STAR_output/*/*bam STAR_bam
	
fi

mkdir STAR_RESULTS
mv STAR_output STAR_bam star_index STAR_RESULTS
#rm -rf STAR-2.5.2b STAR_align_wrapper.sh
#rm -rf "${seqFolder}"
