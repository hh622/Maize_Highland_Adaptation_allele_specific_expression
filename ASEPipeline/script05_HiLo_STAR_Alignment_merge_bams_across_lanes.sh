#!/bin/bash -l
#SBATCH -D /home/hxhu/workdir/scripts
#SBATCH -o /home/hxhu/workdir/logs/05_HiLo_merge_STARrun01_across_lanes/merge_bam_%A_%a.out
#SBATCH -e /home/hxhu/workdir/logs/05_HiLo_merge_STARrun01_across_lanes/merge_bam_%A_%a.err
#SBATCH -J merge_bams
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=48:00:00
#SBATCH -A runciegrp 
#SBATCH -p med2

module load R
module load maven
module load java
module load GATK/4.0
module load samtools

picardDir=/home/hxhu/programs

#reference genome
FASTA=/group/runciegrp/SharedResources/Genomes/Zea_mays/B73/AGPv4.36/Zea_mays.AGPv4.dna.toplevel.fa
DIR_bam_files=/home/hxhu/workdir/04_map2B73v4STAR/pass2NewName
DIR_merge_lanes=/home/hxhu/workdir/05_STARrun01_merge_bams/merge_bams

#get FileNames and SAMPLENAME
FileNames=$(basename -a $DIR_bam_files/*MQ40_rgadded_sorted.bam | cut -d"." -f 1 | cut -f1-6 -d"_" | sort -u)
FILES=(`echo $FileNames`)
echo ${#FILES[@]} #624, length of array

# grab out filename from the array exported from our ‘parent’ shell
SAMPLENAME=${FILES[$SLURM_ARRAY_TASK_ID]}

#-------------------------------
# merge_lanes_filter_dups

echo "step 1: Merge bam files over lanes to sample/library level (1 lib/sample for our data)"
java -jar $picardDir/picard_2.21.7.jar MergeSamFiles $(printf 'I=%s ' $DIR_bam_files/$SAMPLENAME*MQ40_rgadded_sorted.bam) \
        o=$DIR_merge_lanes/$SAMPLENAME.lanemerged.bam

echo "step 2: Filter duplicate reads afer lane merge"
## for paired-end reads: python rmdup_pe.py <sorted.input.bam> <output.bam>
#python $WASP/mapping/rmdup_pe.py $DIR_merge_lanes/$SAMPLENAME.lanemerged.bam $DIR_merge_lanes/$SAMPLENAME.lanemerged.dedupped.bam

echo "step 3: Sort bam"
samtools sort -o $DIR_merge_lanes/$SAMPLENAME.lanemerged.sort.bam \
              $DIR_merge_lanes/$SAMPLENAME.lanemerged.bam

echo "step 4: index sorted dedeupped bam"
samtools index $DIR_merge_lanes/$SAMPLENAME.lanemerged.sort.bam

