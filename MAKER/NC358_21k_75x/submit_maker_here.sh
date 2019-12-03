#!/bin/bash
# specify BASH shell
#$ -S /bin/bash
# run job in the current working directory where qsub is executed from
#$ -cwd
# Tell the SGE that this is an array job, with "tasks" to be numbered 
#  specify that the job requires 8GB of memory for each process
#$ -t 3
#$ -pe mpi 32 -l m_mem_free=4G
# make sure that there is at least 20G of tmp space on the node
#$ -l tmp_free=300G

#module load openmpi-x86_64
source /sonas-hs/ware/hpc/home/kchougul/.bash_profile
export LD_PRELOAD=/sonas-hs/ware/hpc/home/mcampbel/applications/openmpi-1.8.8/build/lib/libmpi.so
export AUGUSTUS_CONFIG_PATH=/sonas-hs/ware/hpc_norepl/data/programs/augustus/augustus.2.5.5/config/

#let i=$SGE_TASK_ID-1
/sonas-hs/ware/hpc/home/mcampbel/applications/openmpi-1.8.8/build/bin/mpiexec -mca btl ^openib -n 32 /sonas-hs/ware/hpc/home/mcampbel/applications/maker/bin/maker -g genome_fasta/contig_fasta/Contig$SGE_TASK_ID.fasta -t 10 -fix_nucleotides -TMP $TMPDIR
