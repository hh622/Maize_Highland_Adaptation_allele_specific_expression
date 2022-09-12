#!/bin/bash -l
#SBATCH -D /home/hxhu/workdir/scripts
#SBATCH -o /home/hxhu/workdir/logs/04_STAR_Alignment_B73v4/staraln_%A_%a.out
#SBATCH -e /home/hxhu/workdir/logs/04_STAR_Alignment_B73v4/staraln_%A_%a.err
#SBATCH -J STAR_alignment
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=240:00:00
#SBATCH -A runciegrp 
#SBATCH -p med2

module load star/2.7.2a

refline=B73
DIR_fastq=/home/hxhu/workdir/01_Fastq
DIR_TrimmomaticQC=/home/hxhu/workdir/03_TrimmomaticQC
DIR_STAR=/home/hxhu/workdir/04_map2B73v4STAR

#preparateion: get bash arry of fastq file names
FileNames=$(basename -a $DIR_fastq/*fq.gz | cut -d"." -f1 | cut -d"_" -f1-8 | sort -u)
FILES=(`echo $FileNames`)
echo ${#FILES[@]} #2354

#grab out each filename from the filenames array(i.e.FILES)
FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]} 
#FILENAME=${FILES[1000]}; echo $FILENAME #e.g. FILENAME=Mex-30-Low_pv_blk2_V4leaftip_batch2_plate2_lane1_ID467

echo "##SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}; FILENAME=${FILENAME}"

#------------------------------------------
#STAR ALIGNMENT
cpunum=10

#reference genome
FASTA=/group/runciegrp/SharedResources/Genomes/Zea_mays/B73/AGPv4.36/Zea_mays.AGPv4.dna.toplevel.fa

DIR_STARpass1=$DIR_STAR/pass1; #mkdir $DIR_STARpass1
DIR_STARpass2=$DIR_STAR/pass2; #mkdir $DIR_STARpass2
genomeDirpass2_all=$DIR_STAR/stargenome_CML228v1.0_pass2; #mkdir $genomeDirpass2_all

echo "step 1: generate STAR genome index for 1-pass alignment (pre-generated)"
genomeDirpass1=/home/hxhu/workdir/genome/STAR_genome_index/${refline}v4

echo "step 2: 1-pass STAR reads alignment"
STAR --genomeDir $genomeDirpass1 --readFilesCommand zcat --readFilesIn $DIR_TrimmomaticQC/${FILENAME}_PE_1.fastq.gz $DIR_TrimmomaticQC/${FILENAME}_PE_2.fastq.gz \
     --outFileNamePrefix $DIR_STARpass1/$FILENAME. --runThreadN $cpunum
     #--outReadsUnmapped $DIR_STARpass1/unmapped.$FILENAME #notes: I don't need to output unmapped reads, so take this option out

echo "step 3: create a new genome index for the 2-pass STAR" #create dir for generate sample-specific genome index for the 2-pass STAR
genomeDirpass2=$genomeDirpass2_all/$FILENAME; mkdir $genomeDirpass2

STAR --runMode genomeGenerate --genomeDir $genomeDirpass2 --genomeFastaFiles $FASTA \
    --sjdbFileChrStartEnd $DIR_STARpass1/$FILENAME.SJ.out.tab --sjdbOverhang 100 --runThreadN $cpunum --limitGenomeGenerateRAM 40000000000
    
echo "step 4: 2-pass STAR reads alignment"
STAR --genomeDir $genomeDirpass2 --readFilesCommand zcat --readFilesIn $DIR_TrimmomaticQC/${FILENAME}_PE_1.fastq.gz $DIR_TrimmomaticQC/${FILENAME}_PE_2.fastq.gz \
    --outFileNamePrefix $DIR_STARpass2/$FILENAME. --outReadsUnmapped Fastx --outSAMmapqUnique 60 --runThreadN $cpunum \
    --outSAMtype BAM SortedByCoordinate

echo "remove STAR genome index for pass2 immediately after the 2nd-pass alignment, which is 20 GB per sample"
rm -r $genomeDirpass2

echo "step 5: Filter_STAR_alignment with MQ>=40"

mv $DIR_STARpass2/$FILENAME.Aligned.sortedByCoord.out.bam $DIR_STARpass2/$FILENAME.sorted.bam
samtools view -b -q 40 $DIR_STARpass2/$FILENAME.sorted.bam > $DIR_STARpass2/${FILENAME}_sorted_MQ40.bam
samtools index $DIR_STARpass2/${FILENAME}_sorted_MQ40.bam

echo "step 6: Addreadgroups and sort bam file according to coordinate"
DIR_Picard=/home/hxhu/programs

java -jar $DIR_Picard/picard_2.21.7.jar AddOrReplaceReadGroups \
        I=$DIR_STARpass2/${FILENAME}_sorted_MQ40.bam \
        O=$DIR_STARpass2/${FILENAME}_MQ40_rgadded_sorted.bam \
        SO=coordinate \
        RGID=$FILENAME \
        RGLB=$(echo $FILENAME | cut -d"_" -f1-4) \
        RGPL=Illumina \
        RGPU=$FILENAME \
        RGSM=$(echo $FILENAME | cut -d"_" -f1-4)

echo "step 7: remove intermediate files/dirs"
rm -r $DIR_STARpass1/$FILENAME*
rm $DIR_STARpass2/$FILENAME.sorted.bam
rm $DIR_STARpass2/${FILENAME}_sorted_MQ40.bam
rm $DIR_STARpass2/${FILENAME}_sorted_MQ40.bam.bai
rm -r $DIR_STARpass2/${FILENAME}._STARtmp
rm $DIR_STARpass2/${FILENAME}.SJ.out.tab
