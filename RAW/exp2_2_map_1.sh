#!/bin/bash

#SBATCH --qos castles
#SBATCH --ntasks 72 # request 8 cores for the job. N.B. check whether the fastq-dump can parallelise, else this is redundant and you should set to "1"
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 3000 # this requests 2 hours, but you will need to adjust depending on runtime. Test job execution time with just a couple of input files then scale accordingly


module purge; module load bluebear
module load bear-apps/2022b
module load ShortStack/4.0.2-foss-2022b

genome="../2_references/SL4_SM41_merged.fa"
genome_index="../2_references/SL4_SM41_merged"

#### STEP B - aligment for de movo sRNA loci identification
#index reference genome
#bowtie-build --threads 6 $genome $genome_index

# this example is for the random format data set
mkdir -p ./3_de_novo_detection
mkdir -p ./3_de_novo_detection/
mkdir -p ./3_de_novo_detection/1_alignment


# map raw files to appropriate genomes
FILES=./2_trimmed/*trimmed.fq

# alignemnt
for f in $FILES
do
f=${f##*/}
f=${f%_trimmed.fq}
#map unique
echo "mapping with bowtie ... $f" >> ./3_de_novo_detection/stats_alignment.txt

(bowtie -v 0 -a -m 1 --sam -x $genome_index "./2_trimmed/$f""_trimmed.fq" | samtools view -bS | samtools sort -o "./3_de_novo_detection/1_alignment/$f"".bam") 2>> ./3_de_novo_detection/stats_alignment.txt

# cluster analysis with ShortStack - use the alignment files as input
echo "clustering with shortstack ... $f" >> ./3_de_novo_detection/stats_alignment.txt

ShortStack \
--bamfile "./3_de_novo_detection/1_alignment/$f"".bam" \
--genomefile $genome \
--threads 6 \
--pad 200 \
--mincov 0.5 \
--nohp \
--outdir ./3_de_novo_detection/$f  >> ./3_de_novo_detection/stats_alignment.txt 2>&1

done


mkdir -p ./4_sRNA_loci
