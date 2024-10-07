#!/bin/bash

## Testing functions 

chmod +x ./PACKAGES/mobileRNAcl/mobileRNAcl.sh
cd ./PACKAGES/mobileRNAcl


## Test help 
./mobileRNAcl.sh RNAmergeGenomes --help

## Test merging genome 
./mobileRNAcl.sh RNAmergeGenomes -i ./tests/Genomes/reduced_chr2_Tomato.fa.gz ./tests/Genomes/reduced_chr12_Eggplant.fa.gz -o ./tests/Genomes/mergedgenometest.fa
./mobileRNAcl.sh RNAmergeAnnotations -i ./tests/Genomes/reduced_chr2_Tomato.gff.gz ./tests/Genomes/reduced_chr12_Eggplant.gff.gz -o ./tests/Genomes/mergedannotationtest.gff



## Test mapping sRNA -- data taken from mobileRNA package
./mobileRNAcl.sh map_sRNA -f ./tests/Genomes/mergedgenometest3.fa -x ./tests/Genomes/mergedgenometest3 -i ./tests/sRNAtestdata/ -o ./tests/sRNAtestoutput/ 


## Test mapping mRNA 

### 1 - make simulated data >> see ./Generating_test_data/
### 2 - convert fasta files to fq 
for f in ./tests/mRNAtestdata/*.fasta; do 
    sample_name=$(basename "$f" | sed -E 's/(\.fasta)$//')
    seqtk seq -F 'I' $f  > ./tests/mRNAtestdata/$sample_name"".fq
done 
rm ./tests/mRNAtestdata/*.fasta




./mobileRNAcl.sh map_mRNA -f ./tests/Genomes/mergedgenometest3.fa -x ./tests/Genomes/mergedgenometest3 -g ./tests/Genomes/mergedgenometest3.gff -i ./tests/mRNAtestdata/ -o ./tests/mRNAtestoutput/ 
