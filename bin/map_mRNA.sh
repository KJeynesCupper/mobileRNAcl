#!/bin/bash

# Function to perform alignment and clustering
map_sRNA() {
    local fasta_dir=""
    local index_dir=""
    local gff_dir=""
    local fastq_dir=""
    local output_dir=""
    local threads=6  
    local paired=FALSE 
    local format=bam 
    local a=0 
    local order=pos 
    local stranded=no 
    local mode=union 
    local nonunique=none 
    local type=mRNA 
    local idattr=Name

    # Part 1 
   
    # Parse command-line options
    while getopts "FASTA:index:GFF:i:o:threads:paired:format:a:order:stranded:mode:nonunique:type:idattr:" opt; do
        case ${opt} in
            FASTA) fasta_dir="$OPTARG" ;;
            index) index_dir="$OPTARG" ;;   
            GFF) gff_dir="$OPTARG" ;;   
            i) fastq_dir="$OPTARG" ;;    
            o) output_dir="$OPTARG" ;;   
            threads) threads="$OPTARG" ;;      
            paired) paired="$OPTARG" ;;
            format) format="$OPTARG" ;;  
            a) a="$OPTARG" ;;  
            order) order="$OPTARG" ;;  
            stranded) stranded="$OPTARG" ;;  
            mode) mode="$OPTARG" ;;  
            nonunique) nonunique="$OPTARG" ;;  
            type) type="$OPTARG" ;;  
            idattr) idattr="$OPTARG" ;;           
            *) echo "Invalid option"; exit 1 ;;
        esac
    done

    # Check if all required arguments are provided
    if [[ -z "$fasta_dir" || -z "$index_dir" || -z "$fastq_dir" || -z "$output_dir" ]]; then
        echo "Usage: map_sRNA -r /path/to/fasta/dir -i /path/to/index/dir -f /path/to/fastq/dir -o /path/to/output/dir [-t threads] [-p pad] [-m mincov] [-a dicermin] [-b dicermax]"
        exit 1
    fi

    # check if hisat, samtools, and htseq are installed: 


    # Check if Bowtie index exists, and build if necessary
    if [[ ! -f "${index_dir}.ht2" | ! -f "${index_dir}.ht2l" ]]; then
        echo "HISAT2 index not found, building index..."
        hisat2-build -p "$threads" "$fasta_dir" "$index_dir"
    fi

    # List all FASTQ files in the directory

    echo "Begining alignment..."
    if [[ "$paired" == "TRUE" ]]; then
    
        local files=$(find "$fastq_dir" -type f \( -name "*_L1_1.fastq" -o -name "*_L1_1.fq" -o -name "*_L1_1.fastq.gz" -o -name "*_L1_1.fq.gz" \))

        for f in $files; do
            local sample_name=$(basename "$f" | sed -E 's/(\_L1_1.fastq|\_L1_1.fq|\_L1_1.fastq\.gz|\_L1_1.fq\.gz)$//')
            local dir_name=$(dirname "$f")
            file1=$dir_name/${sample_name}_L1_1*
            file2=$dir_name/${sample_name}_L1_2*
            file3=$dir_name/${sample_name}_L2_1*
            file4=$dir_name/${sample_name}_L2_2*
            
            mkdir -p "$output_dir/$sample_name"

            hisat2 -p 6 \
            -x $genome_index \
            --summary-file $stats/${f}_summary.txt \
            -1 $dir_name/$file1,$dir_name/$file3 \
            -2 $dir_name/$file2q,$dir_name/$file4 \
            -S $output_dir/$sample_name/${sample_name}_AllReads.sam

            samtools view -h $output_dir/$sample_name/${sample_name}_AllReads.bam > $output_dir/$sample_name/${sample_name}_AllReads.sam
            grep "NH:i:1" $output_dir/$sample_name/${sample_name}_AllReads.sam > $output_dir/$sample_name/${sample_name}_uniqueReads.sam 
            samtools view -H $output_dir/$sample_name/${sample_name}_uniqueReads.sam  > $output_dir/$sample_name/${sample_name}_header.txt
            cat $output_dir/$sample_name/${sample_name}_header.txt $output_dir/$sample_name/${sample_name}_uniqueReads.sam > $output_dir/$sample_name/${sample_name}_uniqueReadsH.sam
            samtools view -Sb $output_dir/$sample_name/${sample_name}_uniqueReadsH.sam | samtools sort -o $output_dir/$sample_name/${sample_name}_unique_sorted.bam
            samtools index $output_dir/$sample_name/${sample_name}_unique_sorted.bam $output_dir/$sample_name/${sample_name}_unique_index.bam
            rm $output_dir/$sample_name/*.sam
            rm $output_dir/$sample_name/*_header.txt

        done
    
    

    else
    
        local files=$(find "$fastq_dir" -type f \( -name "*_1.fastq" -o -name "*_1.fq" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \))

        for f in $files; do
            local sample_name=$(basename "$f" | sed -E 's/(\_1.fastq|\_1.fq|\_1.fastq\.gz|\_1.fq\.gz)$//')
            local dir_name=$(dirname "$f")
            file1=$dir_name/${sample_name}_1*
            file2=$dir_name/${sample_name}_2*
            
            mkdir -p "$output_dir/$sample_name"

            hisat2 -p "$threads" \
            -x $index_dir \
            --summary-file $output_dir/$sample_name/${f}_summary.txt \
            -1 $dir_name/$file1 \
            -2 $dir_name/$file2 \
            -S $output_dir/$sample_name/${sample_name}_AllReads.sam

            samtools view -h $output_dir/$sample_name/${sample_name}_AllReads.bam > $output_dir/$sample_name/${sample_name}_AllReads.sam
            grep "NH:i:1" $output_dir/$sample_name/${sample_name}_AllReads.sam > $output_dir/$sample_name/${sample_name}_uniqueReads.sam 
            samtools view -H $output_dir/$sample_name/${sample_name}_uniqueReads.sam  > $output_dir/$sample_name/${sample_name}_header.txt
            cat $output_dir/$sample_name/${sample_name}_header.txt $output_dir/$sample_name/${sample_name}_uniqueReads.sam > $output_dir/$sample_name/${sample_name}_uniqueReadsH.sam
            samtools view -Sb $output_dir/$sample_name/${sample_name}_uniqueReadsH.sam | samtools sort -o $output_dir/$sample_name/${sample_name}_unique_sorted.bam
            samtools index $output_dir/$sample_name/${sample_name}_unique_sorted.bam $output_dir/$sample_name/${sample_name}_unique_index.bam
            rm $output_dir/$sample_name/*.sam
        done 
    fi 

echo "Begining raw count estimation..." 
local alignmentfiles="$output_dir/$sample_name"

    for i in $alignmentfiles; do
        file=$alignmentfiles/*_unique_sorted.bam
        j=${file##*/}
        sample_name=${j%*}

        python -m HTSeq.scripts.count \
        --format="$format" -a="$a" --order="$order" --stranded="$stranded" --mode="$mode" --nonunique="$nonunique" --type="$type" --idattr="$idattr" \
        $i $GFF > $output_dir/$sample_name/Results.txt
    done 

echo "Complete"

}












