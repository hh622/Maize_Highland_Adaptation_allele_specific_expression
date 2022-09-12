#!/bin/bash -l
#SBATCH -D /home/hxhu/workdir/scripts
#SBATCH -o /home/hxhu/workdir/logs/07_Filter_WGS_VCF/FilterVCF_%A_%a.out
#SBATCH -e /home/hxhu/workdir/logs/07_Filter_WGS_VCF/FilterVCF_%A_%a.err
#SBATCH -J Filter_VCF
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=48:00:00
#SBATCH -A runciegrp 
#SBATCH -p med2

DIR_bam_files=/home/hxhu/workdir/04_map2B73v4STAR/pass2NewName
DIR_reads_counts=/home/hxhu/workdir/06_ASEReadCounter_run01/counts_table
DIR_VCF=/home/hxhu/workdir/01_ANGSD/ANGSD_WGS_LR_wholegenome
DIR_FilteredVCF=/home/hxhu/workdir/07_Filter_WGS_VCF
#mkdir $DIR_FilteredVCF

#get FileNames and SAMPLENAME
FileNames=$(basename -a $DIR_bam_files/*MQ40_rgadded_sorted.bam | cut -d"." -f 1 | cut -f1-6 -d"_" | sort -u)
FILES=(`echo $FileNames`)
echo ${#FILES[@]} #624, length of array

# grab out filename from the array exported from our 'parental' shell
SAMPLENAME=${FILES[$SLURM_ARRAY_TASK_ID]}

#---------------------------------------
#filtering ASEReadCounter results

# filter for 1+REF, 1+ALT, count >= 10 and abs(log2(REF/ALT)) <= 2
tail -n +2 $DIR_reads_counts/${SAMPLENAME}_ASE_read_counts.table | \
	awk '{Thr_totalCount=10; if ($6 >0 && $7 >0 && $8 >= Thr_totalCount) print $0}' OFS='\t' | \
	awk '{ if( log($6/$7)/log(2)<=2 && log($6/$7)/log(2) >= -2 ) print $1,$2 }' OFS='\t' | \
	awk '{$3=$1 "." $2; print $3}' OFS='\t' \
	> $DIR_reads_counts/${SAMPLENAME}_chr_Position_Selected_Thr10_ASE2.txt # wc -l 89458

#subset the raw VCF to make RNAseq sample-specific VCF
awk '{$11=$1 "." $2; print $0}' OFS='\t' $DIR_VCF/WGSLR_wholegenome_phasedGT.vcf |\
       	grep -Fwf $DIR_reads_counts/${SAMPLENAME}_chr_Position_Selected_Thr10_ASE2.txt - | \
	cut -f1-10 | \
	sed '1i #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPhasedGenotype' \
	> $DIR_FilteredVCF/${SAMPLENAME}_WGSLR_wholegenome_phasedGT.vcf #89459 with header

#---------------------------------------------------
#For plotting ASE=refCount/altCount distribution
#awk '{if ($6 >0 && $7 >0) print $0}' OFS='\t' $DIR_reads_counts/${SAMPLENAME}_ASE_read_counts.table \
#	> $DIR_reads_counts2/${SAMPLENAME}_ASE_read_counts.table


