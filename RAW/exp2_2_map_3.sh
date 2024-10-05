#!/bin/bash

#SBATCH --qos castles
#SBATCH --ntasks 72 # request 8 cores for the job. N.B. check whether the fastq-dump can parallelise, else this is redundant and you should set to "1"
#SBATCH --time 5000 # this requests 2 hours, but you will need to adjust depending on runtime. Test job execution time with just a couple of input files then scale accordingly

module purge; module load bluebear
module load bear-apps/2022b
module load ShortStack/4.0.2-foss-2022b


genome="../2_references/SL4_SM41_merged.fa"

mkdir -p ./5a_clustering_EggEgg_vs_EggTom
mkdir -p ./5b_clustering_TomTom_vs_TomEgg


#################################################################################
############ EggEgg_vs_EggTom   ############
#################################################################################

# location of alignment files
FILES="./3_de_novo_detection/1_alignment/EggEgg_1.bam ./3_de_novo_detection/1_alignment/EggEgg_2.bam ./3_de_novo_detection/1_alignment/EggEgg_3.bam ./3_de_novo_detection/1_alignment/EggEgg_4.bam ./3_de_novo_detection/1_alignment/EggTom_1.bam ./3_de_novo_detection/1_alignment/EggTom_2.bam ./3_de_novo_detection/1_alignment/EggTom_3.bam"

# alignemnt
for f in $FILES
do
f=${f##*/}
f=${f%.bam}

echo "clustering with shortstack ... $f" >> ./5a_clustering_EggEgg_vs_EggTom/stats_alignment.txt

ShortStack \
--bamfile "./3_de_novo_detection/1_alignment/$f"".bam" \
--genomefile $genome \
--locifile ./4_sRNA_loci/EggEgg_vs_EggTom_locifile.txt \
--threads 6 \
--pad 200 \
--mincov 0.5 \
--nohp \
--outdir ./5a_clustering_EggEgg_vs_EggTom/$f  >> ./5a_clustering_EggEgg_vs_EggTom/stats_alignment.txt 2>&1

done




#################################################################################
############ TomTom_vs_TomEgg   ############
#################################################################################

# location of alignment files
FILES="./3_de_novo_detection/1_alignment/TomTom_1.bam ./3_de_novo_detection/1_alignment/TomTom_2.bam ./3_de_novo_detection/1_alignment/TomTom_3.bam ./3_de_novo_detection/1_alignment/TomEgg_1.bam ./3_de_novo_detection/1_alignment/TomEgg_2.bam ./3_de_novo_detection/1_alignment/TomEgg_3.bam ./3_de_novo_detection/1_alignment/TomEgg_4.bam"

# alignemnt
for f in $FILES
do
f=${f##*/}
f=${f%.bam}

echo "clustering with shortstack ... $f" >> ./5b_clustering_TomTom_vs_TomEgg/stats_alignment.txt

ShortStack \
--bamfile "./3_de_novo_detection/1_alignment/$f"".bam" \
--genomefile $genome \
--locifile ./4_sRNA_loci/TomTom_vs_TomEgg_locifile.txt \
--threads 6 \
--pad 200 \
--mincov 0.5 \
--nohp \
--outdir ./5b_clustering_TomTom_vs_TomEgg/$f  >> ./5b_clustering_TomTom_vs_TomEgg/stats_alignment.txt 2>&1

done
