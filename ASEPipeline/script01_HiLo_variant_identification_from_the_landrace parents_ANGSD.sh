#!/bin/bash -l
#SBATCH -D /home/hxhu/workdirscripts
#SBATCH -o /home/hxhu/workdirlogs01_ANGSD_WGS_LR_parallel_%A_%a.out
#SBATCH -e /home/hxhu/workdirlogs/01_ANGSD_WGS_LR_parallel_%A_%a.err
#SBATCH -J ANGSD_WGS_LR
#SBATCH -c 20
#SBATCH -t 96:00:00
#SBATCH --mem=40G
#SBATCH -A runciegrp 
#SBATCH -p high

module load angsd

genome=/group/runciegrp/SharedResources/Genomes/Zea_mays/B73/AGPv4.36/Zea_mays.AGPv4.dna.toplevel.fa
DIR_output=/home/hxhu/workdir/ANGSD_WGS_LR_wholegenome

#make bam file list
DIR_WGS_bam=/group/runciegrp/Projects/HiLo_WGS/mapped_reads
ls $DIR_WGS_bam/WGS_LR_P*/*sort.bam | egrep -v "WGS_LR_P2_1_sort.bam|WGS_LR_P2_3_sort.bam" > bam2.filelist

#notes: I excluded the above two fam files becasue they are liekly duplicates, which are from _old folder
#ls $DIR_WGS_bam/WGS_LR_P*/*sort.bam | egrep "WGS_LR_P2_1_sort.bam|WGS_LR_P2_3_sort.bam"
#/group/runciegrp/Projects/HiLo_WGS/mapped_reads/WGS_LR_P2_1_old/WGS_LR_P2_1_sort.bam
#/group/runciegrp/Projects/HiLo_WGS/mapped_reads/WGS_LR_P2_3_old/WGS_LR_P2_3_sort.bam

angsd -GL 2 -P 20 \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 \
-C 50  -minMapQ 20 -minQ 20 -SNP_pval 1e-6 \
-ref $genome  -anc $genome \
-doMaf 2 -doMajorMinor 4 -doSaf 1 \
-r ${SLURM_ARRAY_TASK_ID}: \
-bam bam2.filelist \
-out $DIR_output/WGSLR_chr${SLURM_ARRAY_TASK_ID}

echo "angsd is done!"
#rm  bam2.filelist
