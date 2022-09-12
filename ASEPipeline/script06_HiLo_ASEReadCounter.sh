#!/bin/bash -l
#SBATCH -D /home/hxhu/workdir/scripts
#SBATCH -o /home/hxhu/workdir/logs/06_ASEReadCounter_run01/ASEReadCounter_%A_%a.out
#SBATCH -e /home/hxhu/workdir/logs/06_ASEReadCounter_run01/ASEReadCounter_%A_%a.err
#SBATCH -J ASEReadCounter
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=48:00:00
#SBATCH -A runciegrp 
#SBATCH -p med2

module load R
module load maven
module load java
module load GATK/4.0
module load samtools

#reference genome
FASTA=/group/runciegrp/SharedResources/Genomes/Zea_mays/B73/AGPv4.36/Zea_mays.AGPv4.dna.toplevel.fa
DIR_bam_files=/home/hxhu/workdir/04_map2B73v4STAR/pass2NewName
DIR_merge_lanes=/home/hxhu/workdir/05_STARrun01_merge_bams/merge_bams
DIR_reads_counts=/home/hxhu/workdir/06_ASEReadCounter_run01/counts_table
#mkdir $DIR_merge_lanes
#mkdir $DIR_reads_counts

#get FileNames and SAMPLENAME
FileNames=$(basename -a $DIR_bam_files/*MQ40_rgadded_sorted.bam | cut -d"." -f 1 | cut -f1-6 -d"_" | sort -u)
FILES=(`echo $FileNames`)
echo ${#FILES[@]} #624, length of array

# grab out filename from the array exported from our 'parent' shell
SAMPLENAME=${FILES[$SLURM_ARRAY_TASK_ID]}

#--------------------------------
# reads counts
DIR_VCF=/home/hxhu/workdir/03_ASEPilot/32_ANGSD/ANGSD_WGS_LR_wholegenome
DIR_VCF_meta=/home/hxhu/workdir/03_ASEPilot/32_ANGSD/VCF_meta

#add metainfo to the $VCFfile
#cat $DIR_VCF_meta/B73v4_GATKvarcall_meta_info.vcf $DIR_VCF/WGSLR_wholegenome_phasedGT_exonic.recode.vcf > \
#        $DIR_VCF/WGSLR_wholegenome_phasedGT_with_meta_info.vcf
#notes: already generated it, no need to redo

VCFfile=$DIR_VCF/WGSLR_wholegenome_phasedGT_with_meta_info.vcf

#index vcf
#bgzip -c $VCFfile > ${VCFfile}.gz; tabix -p vcf ${VCFfile}.gz
#notes: already generated it, no need to redo

#reads count with ASEReadCounter
gatk ASEReadCounter -R $FASTA \
        -I $DIR_merge_lanes/${SAMPLENAME}.lanemerged.sort.bam \
        -V ${VCFfile}.gz \
        -DF NotDuplicateReadFilter \
        --output-format TABLE \
        -O $DIR_reads_counts/${SAMPLENAME}_ASE_read_counts.table

#Like most GATK tools, this tools filters out duplicate reads by default. 
#However, some ASE methods recommend including duplicate reads in the analysis, 
#so the DuplicateRead filter can be disabled using the "-DF NotDuplicateReadFilter" flag in the command-line.
