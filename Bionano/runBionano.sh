#!/bin/bash
if [ "$#" -ne 4 ] ; then
echo "please provide:"
echo -e "\t\t(1) scaffolds in fasta format (input assembly)"
echo -e "\t\t(2) bionano map (de novo assembled) in cmap format"
echo -e "\t\t(3) filtered bionano raw molecules in bnx format"
echo -e "\t\t(4) autonoise errbin file from bionano output (usually located in output/contigs/auto_noise/autoNoise1.errbin)"
echo "";
echo "./runBionano.sh <scaffolds.fasta> <bionano.cmap> <raw-molecules.bnx> <autonoise.errbin>" ;
echo "";
exit 0;
fi
# activate environment and load required modules
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate bionano
ml r-devtools/1.12.0-py2-r3.4-47jzlan
ml rm perl
ml rm python
# setup variables
Rscript -e '.libPaths("~/R/x86_64-pc-linux-gnu-library/3.4/")'
EBROOTBIONANOSOLVE="/work/LAS/mhufford-lab/arnstrm/PanAnd/Solve3.4_06042019a"
genome=$1
cmap=$2
raw=$3
errbin=$4
config="$EBROOTBIONANOSOLVE/HybridScaffold/06042019/hybridScaffold_DLE1_config.xml"
outdir=$(pwd)/$(basename ${genome%.*})
mkdir -p ${outdir}
# run hybrid assembly
perl $EBROOTBIONANOSOLVE/HybridScaffold/06042019/hybridScaffold.pl \
-n ${genome} \
-b ${cmap} \
-c ${config} \
-r $EBROOTBIONANOSOLVE/RefAligner/8949.9232rel/sse/RefAligner \
-o ${outdir} \
-f \
-B 2 \
-N 2 \
-y \
-x \
-m ${raw} \
-p $EBROOTBIONANOSOLVE/Pipeline/06042019 \
-q $EBROOTBIONANOSOLVE/RefAligner/8949.9232rel/optArguments_nonhaplotype_noES_noCut_DLE1_saphyr.xml \
-e ${errbin}
