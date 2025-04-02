#!/bin/bash

module load bowtie
module load seqkit
module load pigz
module load samtools
module load trimgalore
module load sratoolkit

export TMPDIR=/lscratch/$SLURM_JOB_ID

SAMPLE=$1 # SRA run ID (e.g., SRR######)
NAME="all_ismA"
IND=/data/jiangak/sir_mapping
THREADS=10 # Number of threads
OUT=/data/jiangak/sir_mapping/ismA_out_prism # path to output directory
TMP=/lscratch/$SLURM_JOB_ID # Path to temporary directory
HUMAN=/gpfs/gsfs12/users/Irp-jiang/share/DB_Share/human/chm13v2.0/chm13v2.0 #bowtie human reference path
BOWTIE2_INDEXES=/data/jiangak/sir_mapping/all_ismA

# Download reads and trim them using trimgalore
echo "Downloading reads"
echo $SAMPLE

# Step 1: Download the SRA file using prefetch
echo "Downloading SRA file using prefetch"
prefetch -v -O $TMP $SAMPLE

# Step 2: Convert the SRA file to FASTQ format using fasterq-dump
echo "Converting SRA to FASTQ using fasterq-dump"
echo fasterq-dump --split-files -e $THREADS -t $TMP -O $TMP ${TMP}/${SAMPLE}/${SAMPLE}.sra*
fasterq-dump --split-files -e $THREADS -t $TMP -O $TMP ${TMP}/${SAMPLE}/${SAMPLE}.sra*


#fasterq-dump -v -e $THREADS -t $TMP -O /lscratch/$SLURM_JOB_ID/ $SAMPLE

echo "Downloaded files:"
ls /lscratch/$SLURM_JOB_ID/

pigz -p $THREADS /lscratch/$SLURM_JOB_ID/${SAMPLE}_1.fastq
pigz -p $THREADS /lscratch/$SLURM_JOB_ID/${SAMPLE}_2.fastq

echo "Trimming reads"
trim_galore -o /lscratch/$SLURM_JOB_ID/ -j $THREADS --paired /lscratch/$SLURM_JOB_ID/${SAMPLE}_1.fastq.gz /lscratch/$SLURM_JOB_ID/${SAMPLE}_2.fastq.gz

#3 Map reads to the human genome reference and remove likely contaminants
bowtie2 -p $THREADS -x $HUMAN -1 /lscratch/$SLURM_JOB_ID/${SAMPLE}_1_val_1.fq.gz \
    -2 /lscratch/$SLURM_JOB_ID/${SAMPLE}_2_val_2.fq.gz \
     | samtools view -bS | samtools fastq -@ $THREADS -f 12 -F 256 \
     -1 /lscratch/$SLURM_JOB_ID/${SAMPLE}_1.fq -2 /lscratch/$SLURM_JOB_ID/${SAMPLE}_2.fq - 

#mv /lscratch/$SLURM_JOB_ID/${SAMPLE}_1_val_1.fq.gz /lscratch/$SLURM_JOB_ID/${SAMPLE}_1.fq.gz
#mv /lscratch/$SLURM_JOB_ID/${SAMPLE}_2_val_2.fq.gz /lscratch/$SLURM_JOB_ID/${SAMPLE}_2.fq.gz

pigz -p $THREADS /lscratch/$SLURM_JOB_ID/${SAMPLE}_1.fq
pigz -p $THREADS /lscratch/$SLURM_JOB_ID/${SAMPLE}_2.fq


# Summarize total number of reads in the QC'd reads
cat /lscratch/$SLURM_JOB_ID/${SAMPLE}_1.fq.gz | seqkit -j $THREADS stats -T > /lscratch/$SLURM_JOB_ID/${SAMPLE}.stats 
cat /lscratch/$SLURM_JOB_ID/${SAMPLE}_2.fq.gz | seqkit -j $THREADS stats -T >> /lscratch/$SLURM_JOB_ID/${SAMPLE}.stats
TOTAL=$(head -n 2 /lscratch/$SLURM_JOB_ID/${SAMPLE}.stats | tail -n 1| cut -f 4)

# Map reads to the bilirubin reductase index
echo "Mapping Reads"
bowtie2 --no-unal -p $THREADS -x ${IND}/${NAME} -1 /lscratch/$SLURM_JOB_ID/${SAMPLE}_1.fq.gz -2 /lscratch/$SLURM_JOB_ID/${SAMPLE}_2.fq.gz \
    | samtools view -bS | samtools sort > /lscratch/$SLURM_JOB_ID/${SAMPLE}.${NAME}.bam
COUNT=$(samtools idxstats /lscratch/$SLURM_JOB_ID/${SAMPLE}.${NAME}.bam | awk '{ SUM += $3} END {print SUM}')

# Print sample, reads mapped to bilirubin reductase, total reads in sample
echo -e "${SAMPLE}\t${COUNT}\t${TOTAL}" > /lscratch/$SLURM_JOB_ID/${SAMPLE}.${NAME}.mapping.summary
mv /lscratch/$SLURM_JOB_ID/${SAMPLE}.${NAME}.bam ${OUT}/
mv /lscratch/$SLURM_JOB_ID/${SAMPLE}.${NAME}.mapping.summary ${OUT}/

rm ${OUT}/${SAMPLE}.${NAME}.bam 
rm /lscratch/$SLURM_JOB_ID/${SAMPLE}*.gz

