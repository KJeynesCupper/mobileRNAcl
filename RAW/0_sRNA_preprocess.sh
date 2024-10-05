#!/bin/bash

# wget https://englandaws-data1.s3-eu-west-1.amazonaws.com/out/CP2022051800020/X204SC23053552-Z01-F003/checkSize.xls
# wget https://englandaws-data1.s3-eu-west-1.amazonaws.com/out/CP2022051800020/X204SC23053552-Z01-F003/MD5.txt
# wget https://englandaws-data1.s3-eu-west-1.amazonaws.com/out/CP2022051800020/X204SC23053552-Z01-F003/X204SC23053552-Z01-F001.tar

# open tar

 tar -xvf *.tar

################################################################
###  Sample files
################################################################
# 1 .after moving raw file, rename

# A511 -
# A512 -
# A513 -
# A514 -
# A515 -
# A516 -
# A517 -
# A518 -
# A519 -
# A520 -
# A521 -
# A522 -
# A523 -
# A524 -


cp ./novogene_files_raw/Exp2_eggplant-tomato-2023-raw-data/sRNAseq_exp2/RAW/X204SC23053552-Z01-F001/01.RawData/A511/*.gz ./Exp2_workplace/sRNAseq/1_raw/TomTom_1.fq.gz
cp ./novogene_files_raw/Exp2_eggplant-tomato-2023-raw-data/sRNAseq_exp2/RAW/X204SC23053552-Z01-F001/01.RawData/A512/*.gz ./Exp2_workplace/sRNAseq/1_raw/TomTom_2.fq.gz
cp ./novogene_files_raw/Exp2_eggplant-tomato-2023-raw-data/sRNAseq_exp2/RAW/X204SC23053552-Z01-F001/01.RawData/A513/*.gz ./Exp2_workplace/sRNAseq/1_raw/TomTom_3.fq.gz
cp ./novogene_files_raw/Exp2_eggplant-tomato-2023-raw-data/sRNAseq_exp2/RAW/X204SC23053552-Z01-F001/01.RawData/A514/*.gz ./Exp2_workplace/sRNAseq/1_raw/TomEgg_1.fq.gz
cp ./novogene_files_raw/Exp2_eggplant-tomato-2023-raw-data/sRNAseq_exp2/RAW/X204SC23053552-Z01-F001/01.RawData/A515/*.gz ./Exp2_workplace/sRNAseq/1_raw/TomEgg_2.fq.gz
cp ./novogene_files_raw/Exp2_eggplant-tomato-2023-raw-data/sRNAseq_exp2/RAW/X204SC23053552-Z01-F001/01.RawData/A516/*.gz ./Exp2_workplace/sRNAseq/1_raw/TomEgg_3.fq.gz
cp ./novogene_files_raw/Exp2_eggplant-tomato-2023-raw-data/sRNAseq_exp2/RAW/X204SC23053552-Z01-F001/01.RawData/A517/*.gz ./Exp2_workplace/sRNAseq/1_raw/TomEgg_4.fq.gz
cp ./novogene_files_raw/Exp2_eggplant-tomato-2023-raw-data/sRNAseq_exp2/RAW/X204SC23053552-Z01-F001/01.RawData/A518/*.gz ./Exp2_workplace/sRNAseq/1_raw/EggEgg_1.fq.gz
cp ./novogene_files_raw/Exp2_eggplant-tomato-2023-raw-data/sRNAseq_exp2/RAW/X204SC23053552-Z01-F001/01.RawData/A519/*.gz ./Exp2_workplace/sRNAseq/1_raw/EggEgg_2.fq.gz
cp ./novogene_files_raw/Exp2_eggplant-tomato-2023-raw-data/sRNAseq_exp2/RAW/X204SC23053552-Z01-F001/01.RawData/A520/*.gz ./Exp2_workplace/sRNAseq/1_raw/EggEgg_3.fq.gz
cp ./novogene_files_raw/Exp2_eggplant-tomato-2023-raw-data/sRNAseq_exp2/RAW/X204SC23053552-Z01-F001/01.RawData/A521/*.gz ./Exp2_workplace/sRNAseq/1_raw/EggEgg_4.fq.gz
cp ./novogene_files_raw/Exp2_eggplant-tomato-2023-raw-data/sRNAseq_exp2/RAW/X204SC23053552-Z01-F001/01.RawData/A522/*.gz ./Exp2_workplace/sRNAseq/1_raw/EggTom_1.fq.gz
cp ./novogene_files_raw/Exp2_eggplant-tomato-2023-raw-data/sRNAseq_exp2/RAW/X204SC23053552-Z01-F001/01.RawData/A523/*.gz ./Exp2_workplace/sRNAseq/1_raw/EggTom_2.fq.gz
cp ./novogene_files_raw/Exp2_eggplant-tomato-2023-raw-data/sRNAseq_exp2/RAW/X204SC23053552-Z01-F001/01.RawData/A524/*.gz ./Exp2_workplace/sRNAseq/1_raw/EggTom_3.fq.gz


# get file integrities to match md5

################################################################
# 2. gunzip
# make sure a "raw" dir is already created containing the OG data
#place raw sample reads into this directory

cd 1_raw
# then unzip and QC
for i in *.gz;
do
gunzip $i;
done

################################################################
# 3. quality control using fastqc
module load FastQC/0.11.9-Java-11
mkdir -p ./qc

for i in *.fq; do fastqc $i -o ./qc/ ; done

 md5sum * > .qc/raw_MD5.txt


################################################################
###    trim raw reads
################################################################
#note: adaper seqeunce needs to be ammended depending on ones used.
source trim.sh


################################################################
###     Mapping 1
################################################################

# pull genome

cp ./workplace-sRNAseq/reference/ref_merged.fasta ./Exp2_workplace/2_references/
source map_1_unique.sh



################################################################
###    Create list of sRNA loci
################################################################
# in R, RNAloci

################################################################
###    Mapping 2
################################################################

source map_2_unique.sh



## done :)
