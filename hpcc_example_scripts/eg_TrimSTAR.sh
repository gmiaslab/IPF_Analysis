#!/bin/bash --login
#
# Example trimming and STAR gene-count generation script.
# This script trims FASTQ files with Trimmomatic and maps reads with STAR using
# --quantMode GeneCounts. STAR gene-count output was then used for downstream
# strandedness/count-column selection.
#
# Input arguments:
# eg_TrimSTAR.sh baseName[1] headcrop[2] crop[3] paired[4] threads[5] cleanup[6] filePattern[7]

####MODULES####
# Example HPCC modules used for the executed workflow. Adjust module names and
# versions for the local cluster environment.
module purge
module load Trimmomatic/0.39-Java-17
module load STAR/2.7.11b-GCC-13.2.0
###*******DEFINITIONS*******###
#working directory
pwd=${pwd}
referenceDirectory="/path/to/reference"
illuminaClipFile="/path/to/TruSeq3-PE.fa"
threads=$5 #number of threads - check that correct no. are allocated in the resources above
filePattern=$7 #pattern for extension
# reference index directory
starIndex="/path/to/STAR_index"
# specify GTF file (for GeneCounts)
gtf="/path/to/gencode.annotation.gtf"
#Select quant mode for either output SAM/BAM or GeneCounts or both
# quantMode="TranscriptomeSAM GeneCounts"
quantMode="GeneCounts"
#Specify kind of output SAM/BAM
SAMType="BAM Unsorted" #(e.g. for SRA results)
# SAMType="BAM SortedByCoordinate" #(e.g. splicing)
# twopassMode="Basic" #(e.g. splicing)
twopassMode="None" #(e.g. for SRA results)

if [ $4 == "yes" ] #paired data
    then
        input1="$1_1${filePattern}"
        input2="$1_2${filePattern}"
        echo "trimming $1" `date`
        # base=$(echo ${input1} | sed 's/.fastq$//g')
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads ${threads} ${input1} ${input2} -baseout "$1.trimmed.fq.gz" ILLUMINACLIP:${referenceDirectory}/TruSeq3-PE.fa:2:30:10 HEADCROP:$2 CROP:$3
        rm $input1
        rm $input2
        # where to write results
        resultsOut="$1_Genes"
        echo "Running Star Paired $1"
        echo `date`
        STAR --runMode alignReads --runThreadN ${threads} --genomeDir ${starIndex} --sjdbGTFfile ${gtf} --quantMode ${quantMode} --outSAMtype ${SAMType} --twopassMode ${twopassMode} --readFilesCommand zcat --outFileNamePrefix ${resultsOut} --readFilesIn $1.trimmed_1P.fq.gz $1.trimmed_2P.fq.gz
        echo "Done Mapping $1"
        echo `date`
fi

if [ $4 == "no" ] #non-paired data
    then
        input1="$1${filePattern}"
        echo "trimming $1" `date`
        # base=$(echo ${input1} | sed 's/.fastq$//g')
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads ${threads} ${input1} "$1.trimmed.fq.gz" ILLUMINACLIP:${referenceDirectory}/TruSeq3-SE.fa:2:30:10 HEADCROP:$2 CROP:$3
        rm $input1
        # where to write results
        resultsOut="$1_Genes"
        echo "Running Star SE $1"
        echo `date`
        STAR --runMode alignReads --runThreadN ${threads} --genomeDir ${starIndex} --sjdbGTFfile ${gtf} --quantMode ${quantMode} --outSAMtype ${SAMType} --twopassMode ${twopassMode} --readFilesCommand zcat --outFileNamePrefix ${resultsOut} --readFilesIn $1.trimmed.fq.gz
echo "Done Mapping $1"
        echo `date`
fi

####*****FILE CLEANUP REMOVAL *****#####
#Comment out if do not want to remove sam and Gene directory
if [ "$6" == "yes" ]
    then
        echo "File Cleanup"
        rm -R "${resultsOut}_STARgenome"
        rm "${resultsOut}Aligned.out.bam"
fi

###*********NOTES*********###
# Trimmomatic ILLUMINACLIP removes adapter and other Illumina-specific sequence.
# HEADCROP removes bases from the 5' end of each read, and CROP truncates reads
# to the specified mapped length. STAR --quantMode GeneCounts produces the
# gene-count table used by eg_strandedness_count_selection.R.
