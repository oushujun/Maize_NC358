module load BLAST/2.2.26-Linux_x86_64
formatdb -p F -i $genome
blastall -p blastn -d $genome -i knob180.fa -m 8 -b 5000 > knob180.blast
blastall -p blastn -d $genome -i TR-1.fa -m 8 -b 5000 > TR-1.blast
blastall -p blastn -d $genome -i CentC.fa -m 8 -b 5000 > CentC.blast
blastall -p blastn -d $genome -i rDNA_intergenic_spacer.fa -m 8 -b 5000 > rDNA_intergenic_spacer.blast
