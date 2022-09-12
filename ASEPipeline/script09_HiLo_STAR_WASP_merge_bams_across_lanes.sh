#!/bin/bash -l
#SBATCH -D /home/hxhu/workdir/12_HiLoASE_B73v4/scripts
#SBATCH -o /home/hxhu/workdir/12_HiLoASE_B73v4/logs/10_STAR-WASP_merge_bams/mergebam_%A_%a.out
#SBATCH -e /home/hxhu/workdir/12_HiLoASE_B73v4/logs/10_STAR-WASP_merge_bams/mergebam_%A_%a.err
#SBATCH -J mergeBAM
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

refgen=B73v4wholegenome

#reference genome
FASTA=/group/runciegrp/SharedResources/Genomes/Zea_mays/B73/AGPv4.36/Zea_mays.AGPv4.dna.toplevel.fa

DIR_STAR_wasp=/home/hxhu/workdir/12_HiLoASE_B73v4/09_STAR_WASP/pass2
DIR_merge_lanes=/home/hxhu/workdir/12_HiLoASE_B73v4/10_STAR-WASP_merge_bams/merged_bams
DIR_readsCount=/home/hxhu/workdir/12_HiLoASE_B73v4/10_STAR-WASP_merge_bams/readsCount
#mkdir $DIR_merge_lanes
#mkdir $DIR_reads_counts

#get FileNames and SAMPLENAME
FileNames=$(basename -a $DIR_STAR_wasp/*_WASPfilterPASS_RG_sort.bam | cut -d"." -f 1 | cut -f1-6 -d"_" | sort -u)
FILES=(`echo $FileNames`)
echo ${#FILES[@]} #624, length of array

# grab out filename from the array exported from our 'parental' shell
SAMPLENAME=${FILES[$SLURM_ARRAY_TASK_ID]}

##-------------------------------
#step 6: merge_lanes_filter_dups

echo "step 6-1: Merge bam files over lanes to sample/library level (1 lib/sample for our data)"
java -jar $picardDir/picard_2.21.7.jar MergeSamFiles $(printf 'I=%s ' $DIR_STAR_wasp/$SAMPLENAME*_WASPfilterPASS_RG_sort.bam) \
        o=$DIR_merge_lanes/$SAMPLENAME.lanemerged.bam

#echo "step 6-2: Filter duplicate reads afer lane merge"
## for paired-end reads: python rmdup_pe.py <sorted.input.bam> <output.bam>
#python $WASP/mapping/rmdup_pe.py $DIR_merge_lanes/$SAMPLENAME.lanemerged.bam $DIR_merge_lanes/$SAMPLENAME.lanemerged.dedupped.bam

echo "step 6-3: Sort bam"
samtools sort -o $DIR_merge_lanes/$SAMPLENAME.lanemerged.sort.bam \
              $DIR_merge_lanes/$SAMPLENAME.lanemerged.bam

echo "step 6-4: index sorted dedeupped bam"
samtools index $DIR_merge_lanes/$SAMPLENAME.lanemerged.sort.bam

echo "delete unsorted bams to save space"
rm $DIR_merge_lanes/$SAMPLENAME.lanemerged.bam

#-------------------------------------------------------
#partition filtered bam to paternal and maternal origins
samtools view -h $DIR_merge_lanes/$SAMPLENAME.lanemerged.sort.bam | egrep "vA:B:c,1|@" | awk '{if ($1~"@" || $12 !~ /,2|,3/)  print $0}' | samtools view -Shb - > \
        $DIR_merge_lanes/$SAMPLENAME.lanemerged.sort.vA1.bam
samtools view -h $DIR_merge_lanes/$SAMPLENAME.lanemerged.sort.bam | egrep "vA:B:c,2|@" | awk '{if ($1~"@" || $12 !~ /,1|,3/)  print $0}' | samtools view -Shb - >  \
        $DIR_merge_lanes/$SAMPLENAME.lanemerged.sort.vA2.bam

#-------------------------------
#readsCount
readsCount_STARwasp_all0=$(samtools view -c $DIR_merge_lanes/$SAMPLENAME.lanemerged.bam)
readsCount_STARwasp_all=$(samtools view -c $DIR_merge_lanes/$SAMPLENAME.lanemerged.sort.bam)
readsCount_STARwasp_vA1=$(samtools view -c $DIR_merge_lanes/$SAMPLENAME.lanemerged.sort.vA1.bam)
readsCount_STARwasp_vA2=$(samtools view -c $DIR_merge_lanes/$SAMPLENAME.lanemerged.sort.vA2.bam)
echo "$SLURM_ARRAY_TASK_ID,$SAMPLENAME,$readsCount_STARwasp_all0,$readsCount_STARwasp_all,$readsCount_STARwasp_vA1,$readsCount_STARwasp_vA2" > \
	$DIR_readsCount/${SAMPLENAME}_readsCount.csv
