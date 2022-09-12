#!/bin/bash -l
#SBATCH -D /home/hxhu/workdir/scripts
#SBATCH -o /home/hxhu/workdir/logs/10_featureCounts/log_%A.out
#SBATCH -e /home/hxhu/workdir/logs/10_featureCounts/log_%A.err
#SBATCH -J featureCounts
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G
#SBATCH --time=48:00:00
#SBATCH -A runciegrp 
#SBATCH -p med2

DIR_B73v4=/group/runciegrp/SharedResources/Genomes/Zea_mays/B73/AGPv4.36
DIR_merge_lanes=/home/hxhu/workdir/09_STAR-WASP_merge_bams/merged_bams

DIR_subread=/home/hxhu/programs/subread-2.0.1/bin
DIR_reads_counts=/home/hxhu/workdir/10_featureCounts

##run featureCounts
$DIR_subread/featureCounts -T 10 -s 1 \
 -p \
 -g gene_id \
 -t exon \
 -a $DIR_B73v4/Zea_mays.AGPv4.36.IsoformswitchanalyzeR.gtf \
 -o $DIR_reads_counts/HiLo_featurecounts_B73v4wholegenome.results \
 $DIR_merge_lanes/*.lanemerged.sort.vA*.bam

#extract read count matrix from original output of featureCounts
cut -f1,7- $DIR_reads_counts/HiLo_featurecounts_B73v4wholegenome.results | \
	grep -v "#" - | sed 's/\/home\/hxhu\/workdir\/10_HiLoASE_B73v4\/09_STAR-WASP_merge_bams\/merged_bams\///g'| \
	sed 's/lanemerged.sort.//g' | \
	sed 's/\.bam//g' \
        > $DIR_reads_counts/HiLo_featurecounts_ALL_B73v4wholegenome_Thr10_ASE2.table
