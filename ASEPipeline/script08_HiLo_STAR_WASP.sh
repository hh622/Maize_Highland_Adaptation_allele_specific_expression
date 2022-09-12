#!/bin/bash -l
#SBATCH -D /home/hxhu/workdir/scripts
#SBATCH -o /home/hxhu/workdir/logs/08_STAR_WASP/STAR_WASP_%A_%a.out
#SBATCH -e /home/hxhu/workdir/logs/08_STAR_WASP/STAR_WASP_%A_%a.err
#SBATCH -J STAR_WASP_B73v4
#SBATCH -c 10
#SBATCH -t 96:00:00
#SBATCH --mem=40G
#SBATCH -A runciegrp 
#SBATCH -p med2

#export PATH=/programs/STAR-2.7.0f/bin/Linux_x86_64:$PATH
module load star
module load samtools

refline=B73
DIR_fastq=/home/hxhu/workdir/01_FastqNewName
DIR_TrimmomaticQC=/home/hxhu/workdir/03_TrimmomaticQCNewName
DIR_STAR=/home/hxhu/workdir/08_STAR_WASP
DIR_FilteredVCF=/home/hxhu/workdir/07_Filter_WGS_VCF
DIR_readsCount=/home/hxhu/workdir/08_STAR_WASP/readsCount
#mkdir $DIR_output
cpunum=10

#preparateion: get bash arry of fastq file names
FileNames=$(basename -a $DIR_fastq/*fq.gz | cut -d"." -f1 | cut -d"_" -f1-9 | sort -u)
FILES=(`echo $FileNames`)
echo ${#FILES[@]} #2354

#grab out each filename from the filenames array(i.e.FILES)
FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]}

#extract SAMPLENAME from FILENAME
SAMPLENAME=$(echo $FILENAME | cut -f1-6 -d"_")

#------------------------------------------
#step 5: STAR ALIGNMENT
cpunum=10

echo "TASK_ID=$SLURM_ARRAY_TASK_ID; FILENAME=$FILENAME; SAMPLENAME=$SAMPLENAME"

#reference genome
FASTA=/group/runciegrp/SharedResources/Genomes/Zea_mays/B73/AGPv4.36/Zea_mays.AGPv4.dna.toplevel.fa

DIR_STARpass1=$DIR_STAR/pass1; #mkdir $DIR_STARpass1
DIR_STARpass2=$DIR_STAR/pass2; #mkdir $DIR_STARpass2
genomeDirpass2_all=$DIR_STAR/stargenome_pass2; #mkdir $genomeDirpass2_all

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

echo "step 4: STAR-WASP (2-pass STAR reads alignment)"

genomeDirpass2=$genomeDirpass2_all/$FILENAME #sample-specific genome index for the 2-pass STAR
VCFfile=$DIR_FilteredVCF/${SAMPLENAME}_WGSLR_wholegenome_phasedGT.vcf
STAR --genomeDir $genomeDirpass2 --readFilesCommand zcat --readFilesIn $DIR_TrimmomaticQC/${FILENAME}_PE_1.fastq.gz $DIR_TrimmomaticQC/${FILENAME}_PE_2.fastq.gz \
    --outFileNamePrefix $DIR_STARpass2/$FILENAME. --outSAMmapqUnique 60 --runThreadN $cpunum \
    --waspOutputMode SAMtag --varVCFfile $VCFfile \
    --outSAMattributes NH HI AS nM vA vG vW \
    --outSAMtype BAM SortedByCoordinate
    #don't need to output Unmapped reads, so take out "--outReadsUnmapped Fastx"

#echo "remove STAR genome index for pass2 immediately after the 2nd-pass alignment, which is 20 GB per sample"
#rm -r $genomeDirpass2

echo "filter bam file for reads passed WASP filtering"
DIR_STAR_wasp=$DIR_STARpass2
samtools view -h $DIR_STAR_wasp/${FILENAME}.Aligned.sortedByCoord.out.bam | \
	awk '{if($1~"@" || $18=="vW:i:1") print $0}' | samtools view -Shb - \
	> $DIR_STAR_wasp/${FILENAME}_WASPfilterPASS.bam

echo "Addreadgroups and sort bam file according to coordinate"
picardDir=/home/hxhu/programs
java -jar $picardDir/picard_2.21.7.jar AddOrReplaceReadGroups \
        I=$DIR_STAR_wasp/${FILENAME}_WASPfilterPASS.bam \
        O=$DIR_STAR_wasp/${FILENAME}_WASPfilterPASS_RG_sort.bam \
        SO=coordinate \
        RGID=$FILENAME \
        RGLB=$(echo $FILENAME | cut -d"_" -f1-6) \
        RGPL=Illumina \
        RGPU=$FILENAME \
        RGSM=$(echo $FILENAME | cut -d"_" -f1-6)

echo "remove un-used intermediate files/dirs"
rm -r $DIR_STARpass1/$FILENAME*
rm $DIR_STAR_wasp/${FILENAME}_WASPfilterPASS.bam
rm -r $DIR_STARpass2/${FILENAME}._STARtmp
rm $DIR_STARpass2/${FILENAME}.SJ.out.tab

echo "readsCount"
DIR_STAR_wasp=/home/hxhu/workdir/08_STAR_WASP/pass2
readsCount_STARwasp_Fastq=$(samtools view -c $DIR_STAR_wasp/${FILENAME}_WASPfilterPASS_RG_sort.bam)
echo "$SLURM_ARRAY_TASK_ID,$FILENAME,$SAMPLENAME,$readsCount_STARwasp_Fastq" > $DIR_readsCount/${FILENAME}_readsCount.csv
