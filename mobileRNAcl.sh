#!/bin/bash

# Source function scripts  
source ./bin/RNAmergeGenomes.sh   
source ./bin/RNAmergeAnnotations.sh
source ./bin/map_sRNA.sh   
source ./bin/map_mRNA.sh   
source ./bin/metrics_sRNA.sh
source ./bin/metrics_mRNA.sh   


#  usage function
usage() {
    case "$1" in
        RNAmergeGenomes)
            echo "Usage: mobileRNAcl RNAmergeGenomes [options]"
            echo "Options:"
            echo "  --help Display this help message for RNAmergeGenomes"
            echo "  -i     the genome reference files seperated by tab (eg. reference1.fa reference2.fa)"
            echo "  -o     the output FASTA file "
            
            echo "Description:"
            echo "  Merges any number of FASTA files into a merged genome reference assembly.
                    It adds a unique identifier to each FASTA file in alphabetical order. Here, the first 
                    FASTA file will have the prefix `A_` added to each chromosome, the secon FASTA will have the 
                    prefix `B_` added to each chromosome, and so on."
            ;;
        RNAmergeAnnotations)
            echo "Usage: mobileRNAcl RNAmergeAnnotations [options]"
            echo "Options:"
            echo "  --help Display this help message for RNAmergeAnnotations"
            echo "  -i     the genome annotation files seperated by tab (eg. reference1.fa reference2.fa)"
            echo "  -o     the output GFF file"
            
            echo "Description:"
            echo "  Merges any number of GFF files into a merged genome reference assembly.
                    It adds a unique identifier to each GFF file in alphabetical order. Here, the first 
                    GFF file will have the prefix `A_` added to each chromosome, the secon GFF will have the 
                    prefix `B_` added to each chromosome, and so on."
            ;;
        map_sRNA)
            echo "Usage: mobileRNAcl map_sRNA [options]"
            echo "Options:"
            echo "  --help          Display this help message for map_sRNA"
            echo "  -FASTA          merged genome reference (FASTA), path to genome reference"
            echo "  -index          bowtie genome reference index, path to index location" 
            echo "  -i              input files (FASTQ), path to sequencing reads" 
            echo "  -o              output location of results, directory to store output."
            echo "  -threads        threads, set the number of threads to use where more threads means a faster completion time (default = 6)"
            echo "  -pad            pad, initial peaks are merged if they are this distance or less from each other. Must >= 1 (default = 200)" 
            echo "  -mincov         mincov, minimum alignment depth, in units of reads per million, required to nucleate a small RNA cluster during de novo cluster search. Must be a number > 0. (default = 0.5)"
            echo "  -dicermin       dicermin, the minimum size in nucleotides of a valid small RNA (default = 20)"
            echo "  -dicermax       dicermax, the maximum size in nucleotides of a valid small RNA (default = 24)"
            echo "  -dn_mirna       dn_mirna, activates a de novo comprehensive genome-wide search for miRNA loci (defalt = FALSE)"
            
            echo "Description:"
            echo "  Aligns small RNA sequencing reads to given genome reference using Bowtie, retaining only uniquely aligned reads. 
                    To detect small RNAs, ShortStack is utilised. The alignment files for each sample are supplied to ShortStack and 
                    the inital de novo assessment is undertaken. For each sample, a GFF file is generated that stores the detected 
                    sRNA-producing genes. These are merged into a single locifile.txt which is utilised in the next step. The final 
                    sRNA clustering analysis is undertaken with ShortStack and the locifile.txt for each sample."
            ;;
        map_mRNA)
            echo "Usage: mobileRNAcl map_sRNA [options]"
            echo "Options:"
            echo "  --help          Display this help message for map_mRNA"
            echo "  -FASTA          merged genome reference (FASTA)"
            echo "  -index          bowtie genome reference index" 
            echo "  -GFF            merged genome annotation (GFF)" 
            echo "  -i              input files (FASTQ), path to sequencing reads" 
            echo "  -o              output location of results"
            echo "  -threads        threads (default = 6)"
            echo "  -paired         paired (default = FALSE)" 
            echo "  -format         format (default = bam)"
            echo "  -a              a (default = 0)"
            echo "  -order          order (default = pos)"
            echo "  -stranded       stranded (default = NO)"
            echo "  -mode           mode (default = union)"
            echo "  -nonunique      nonunique (default = non)"
            echo "  -type           type (default = mRNA)"
            echo "  -idattr         idattr (default = Name)"
            
            echo "Description:"
            echo "  Aligns messenger RNA sequencing reads to given genome reference using HISAT2, retaining only uniquely aligned reads. 
                    Raw count estimation is undertaken using HTseq. "
            ;;
        *)
            echo "Usage: mobileRNAcl [function] [options]"
            echo "function:"
            echo "  RNAmergeGenomes         merge genome references"
            echo "  RNAmergeAnnotations     merge genome annoations"
            echo "  map_sRNA                Map and cluster small RNA sequencing reads"
            echo "  map_mRNA                Map and cluster messenger RNA sequencing reads"
            echo "  metrics_sRNA            csv file of mapping metris for small RNA"
            echo "  metrics_mRNA            csv file of mapping metris for messenger RNA"
            echo "  --help          Display this help message"
            ;;
    esac
    exit 1
}




# Check if no arguments are passed, show general usage
if [ $# -eq 0 ]; then
    usage
fi

# Parse the first argument and dynamically call the corresponding function or display help
case "$1" in
    RNAmergeGenomes)
        # Check for --help flag
        if [ "$2" == "--help" ]; then
            usage RNAmergeGenomes
        else
            echo "Running RNAmergeGenomes..."
            RNAmergeGenomes 
        fi
        ;;
    RNAmergeAnnotations)
        if [ "$2" == "--help" ]; then
            usage RNAmergeAnnotations
        else
            echo "Running RNAmergeAnnotations..."
            RNAmergeAnnotations
        fi
        ;;
    map_sRNA)
        if [ "$2" == "--help" ]; then
            usage map_sRNA
        else
            echo "Running map_sRNA..."
            map_sRNA
        fi
        ;;
    map_mRNA)
        if [ "$2" == "--help" ]; then
            usage map_mRNA
        else
            echo "Running map_mRNA..."
            map_mRNA
        fi
        ;;
    metric_sRNA)
        if [ "$2" == "--help" ]; then
            usage metric_sRNA
        else
            echo "Running metric_sRNA..."
            metric_sRNA
        fi
        ;;
    metric_mRNA)
        if [ "$2" == "--help" ]; then
            usage metric_mRNA
        else
            echo "Running metric_mRNA..."
            metric_mRNA
        fi
        ;;
    --help)
        usage
        ;;
    *)
        echo "Invalid option: $1"
        usage
        ;;
esac
