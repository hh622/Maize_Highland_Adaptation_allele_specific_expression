#!/bin/bash -l
#SBATCH -D /home/hxhu/workdir/scripts
#SBATCH -o /home/hxhu/workdir/logs/03_HiLo_TrimmomaticQC_%A_%a.out
#SBATCH -e /home/hxhu/workdir/logs/03_HiLo_TrimmomaticQC_%A_%a.err
#SBATCH -J script03_HiLo_TrimmomaticQC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=48:00:00
#SBATCH -A runciegrp 
#SBATCH -p med2

module load star/2.7.2a

cml=CML103
DIR_fastq=/home/hxhu/workdir/11_HiLoASE_CMLs/01_Fastq
DIR_TrimmomaticQC=/home/hxhu/workdir/11_HiLoASE_CMLs/03_TrimmomaticQC
DIR_STAR=/home/hxhu/workdir/11_HiLoASE_CMLs/05_map2CML103STAR

#------------------------------------------------
#Step 1: preparateion: get bash arry of fastq file names
FileNames=$(basename -a $DIR_fastq/*fq.gz | cut -d"." -f1 | cut -d"_" -f1-8 | sort -u)
FILES=(`echo $FileNames`)
#echo ${#FILES[@]} #2354
#ls $DIR_fastq/*fq.gz | cut -d"/" -f7 | cut -d"." -f1 | cut -f8 -d"_" | sort -u |wc -l #2354 (ID counts)

#grab out each filename from the filenames array(i.e.FILES)
FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]} 
#FILENAME=${FILES[1000]}; echo $FILENAME #e.g. FILENAME=Mex-30-Low_pv_blk2_V4leaftip_batch2_plate2_lane1_ID467

echo "##SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}; FILENAME=${FILENAME}"

#----------------------------
#Step 2: fastqc before QC
echo "step 2: run fastqc before QC- skip it temperarily......"

#----------------------------
#Step 3: QC with trimmomatic

echo "step 3: QC with trimmomatic through all fastq files"
java -jar ~/programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 $DIR_fastq/${FILENAME}_1.fq.gz  $DIR_fastq/${FILENAME}_2.fq.gz \
$DIR_TrimmomaticQC/${FILENAME}_PE_1.fastq $DIR_TrimmomaticQC/${FILENAME}_SE_1.fastq $DIR_TrimmomaticQC/${FILENAME}_PE_2.fastq $DIR_TrimmomaticQC/${FILENAME}_SE_2.fastq \
ILLUMINACLIP:/group/runciegrp/SharedResources/Adapters/illumina_adapters_bradseq.fa:2:30:10:2:true \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

rm $DIR_TrimmomaticQC/${FILENAME}_SE_1.fastq $DIR_TrimmomaticQC/${FILENAME}_SE_2.fastq


