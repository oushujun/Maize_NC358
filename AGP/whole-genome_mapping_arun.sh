### use the `$target` that is already set
### don't run the mummer part, it takes a long time and it adds nothing to minimap2 alignment
### before this, be sure to run the repeatmasking on your genomes (before aligning them)
### it takes about 2-3hrs for repeatmasking, few mins for minimap2
### and then use dotplotly to make the dotplots (I've scripts for that too in the same dir)
### this is just to make sure that all the scafs are in order (by comparing it to the NAM NC358 genome)
### Arun Seetharam (05/27/2019)

#!/bin/bash
query=$1
target="/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/b-repeatmasking/B73/B73.pseudomolecules-v1.fasta.masked"
bquery=$(basename $query |sed 's/.masked-pseudomolecules.fasta//g')
btarget=$(basename ${target%.*})
module load mummer
out="${bquery}_${btarget}-v5_mummer"
nucmer --maxmatch -l 80 -c 100 ${target} ${query} --prefix=${out}.nucmer
delta-filter -r ${out}.nucmer.delta > ${out}.nucmer.delta.filter
show-coords -c ${out}.nucmer.delta.filter > ${out}.nucmer.delta.filter.coords
