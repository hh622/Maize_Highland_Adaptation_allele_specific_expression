#!/bin/bash -l
#SBATCH -D /home/hxhu/workdir/scripts
#SBATCH -o /home/hxhu/workdir/logs/script32c_generate_STAR-WASP_input_VCF_whole_genome_%A.out
#SBATCH -e /home/hxhu/workdir/logs/script32c_generate_STAR-WASP_input_VCF_whole_genome_%A.err
#SBATCH -J generate_STAR-WASP_input_VCF_whole_genome
#SBATCH -c 2
#SBATCH -t 96:00:00
#SBATCH --mem=10G
#SBATCH -A runciegrp 
#SBATCH -p med2


DIR=/home/hxhu/workdir/01_ANGSD_WGS_LR_wholegenome
cd $DIR

#generate the dots for missing columns/fields
cut -f1 WGSLR_Wholegenome.mafs | sed 's/chromo/EmptyColumn/g' | sed -e 's/[0-9]/./g' | sed 's/\.././g' > WGSLR_Wholegenome_dots.txt
#wc -l WGSLR_Wholegenome_dots.txt #53891496

#generate Field09_Format
sed 's/EmptyColumn/FORMAT/g' WGSLR_Wholegenome_dots.txt | sed 's/\./GT/g' > WGSLR_Wholegenome_Field09_Format.txt

#generate Field10_GT
sed 's/EmptyColumn/PhasedGenotype/g' WGSLR_Wholegenome_dots.txt | sed 's/\./0|1/g' > WGSLR_Wholegenome_Field10_phasedGT.txt
#wc -l WGSLR_Wholegenome_Field10_phasedGT.txt #53891496

#paste/re-organize Fields 1-10 (STAR-WASP requires each VCF has only one sample in column 10, columns 1-9 are meta info)
cut -f1-4 WGSLR_Wholegenome.mafs | paste - WGSLR_Wholegenome_dots.txt | awk '{$6=$5; $7=$5; $8=$5; $5=$4; $4=$3; $3=$6; print $0 }' OFS='\t' | paste - WGSLR_Wholegenome_Field09_Format.txt WGSLR_Wholegenome_Field10_phasedGT.txt | tail -n +2 | sed '1i #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPhasedGenotype' > WGSLR_Wholegenome_phasedGT.vcf
#wc -l WGSLR_Wholegenome_phasedGT.vcf #53891496

#---------------------------------------------
#subset ANGSD_wholegenome VCF for exonic SNPs

#step 1: make exon bed file
module load bedtools/2.27.1
export PATH=/home/hxhu/programs/bedops_v2.4.37/bin:$PATH
DIR_B73v4=/group/runciegrp/SharedResources/Genomes/Zea_mays/B73/AGPv4.36
DIR_B73v4homedir=/home/hxhu/workdir/genome/B73AGPv4

grep "exon" $DIR_B73v4/Zea_mays.AGPv4.36.gff3 | gff2bed | cut -f1-3 | bedtools merge -i - | grep '^[0-9]' | sort -k1,1n -k2,2n > $DIR_B73v4homedir/Zea_mays.AGPv4.36.exons.bed
#cut -f1 $DIR_B73v4homedir/Zea_mays.AGPv4.36.exons.bed | sort -nu #1 2 3 4 5 6 7 8 9 10

#step 2: extract exonic variants
module load vcftools
bedfile=/home/hxhu/workdir/genome/B73AGPv4/Zea_mays.AGPv4.36.exons.bed
vcftools --vcf WGSLR_Wholegenome_phasedGT.vcf --bed $bedfile \
        --out WGSLR_Wholegenome_phasedGT_exonic --recode --keep-INFO-all
