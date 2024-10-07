#!/bin/bash

# Function to perform alignment and clustering
map_mRNA() {
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
    local python=python

    # Part 1 
    while getopts ":f:x:g:i:o:t:p:f:r:s:m:n:j:d:y:" opt; do
            case ${opt} in
                f) fasta_dir="$OPTARG" ;;
                x) index_dir="$OPTARG" ;;
                g) GFF="$OPTARG" ;;
                i) fastq_dir="$OPTARG" ;;
                o) output_dir="$OPTARG" ;;
                t) threads="$OPTARG" ;;
                p) paired="$OPTARG" ;;
                f) format="$OPTARG" ;;
                r) order="$OPTARG" ;;
                s) stranded="$OPTARG" ;;
                m) mode="$OPTARG" ;;
                n) nonunique="$OPTARG" ;;
                j) type="$OPTARG" ;;
                d) idattr="$OPTARG" ;;
                y) python="$OPTARG" ;;
                \?) 
                    echo "Invalid option: -$OPTARG" >&2
                    
                    exit 1
                    ;;
                :) 
                    echo "Option -$OPTARG requires an argument." >&2
                    
                    exit 1
                    ;;
            esac
        done


echo ""
echo "FASTA: $fasta_dir"
echo "GFF: $GFF"
echo "Input: $fastq_dir"
echo "Output: $output_dir"
echo "Index: $index_dir"
echo "Threads: $threads"
echo "Paired: $paired"
echo "Format: $format"
echo "Order: $order"
echo "Stranded: $stranded"
echo "Mode: $mode"
echo "Nonunique: $nonunique"
echo "Type: $type"
echo "idattr: $idattr"
echo "python: $python"
echo ""



    echo "HISAT2 index not found, building index..."
    hisat2-build -p "$threads" "$fasta_dir" "$index_dir"


    # List all FASTQ files in the directory

    echo "Begining alignment..."
    if [[ "$paired" == "TRUE" ]]; then
    
        local files=$(find "$fastq_dir" -type f \( -name "*_L1_1.fastq" -o -name "*_L1_1.fq" -o -name "*_L1_1.fastq.gz" -o -name "*_L1_1.fq.gz" \))

        for f in $files; do
            local sample_name=$(basename "$f" | sed -E 's/(\_L1_1.fastq|\_L1_1.fq|\_L1_1.fastq\.gz|\_L1_1.fq\.gz)$//')
            local dir_name=$(dirname "$f")

            if [[ $f == *".fq" ]]; then
            file1=$dir_name/${sample_name}_L1_1.fq
            file2=$dir_name/${sample_name}_L1_2.fq
            file3=$dir_name/${sample_name}_L2_1.fq
            file4=$dir_name/${sample_name}_L2_2.fq
            fi
            
            if [[ $f == *".fq.gz" ]]; then
            file1=$dir_name/${sample_name}_L1_1.fq.gz
            file2=$dir_name/${sample_name}_L1_2.fq.gz
            file3=$dir_name/${sample_name}_L2_1.fq.gz
            file4=$dir_name/${sample_name}_L2_2.fq.gz
            fi

            if [[ $f == *".fastq" ]]; then
            file1=$dir_name/${sample_name}_L1_1.fastq
            file2=$dir_name/${sample_name}_L1_2.fastq
            file3=$dir_name/${sample_name}_L2_1.fastq
            file4=$dir_name/${sample_name}_L2_2.fastq
            fi
            
            if [[ $f == *".fastq.gz" ]]; then
            file1=$dir_name/${sample_name}_L1_1.fastq.gz
            file2=$dir_name/${sample_name}_L1_2.fastq.gz
            file3=$dir_name/${sample_name}_L2_1.fastq.gz
            file4=$dir_name/${sample_name}_L2_2.fastq.gz
            fi

            
            
            mkdir -p "$output_dir/$sample_name"

            hisat2 -p 6 \
            -x $genome_index \
            --summary-file $stats/${f}_summary.txt \
            -1 $file1,$file3 \
            -2 $file2,$file4 \
            -S $output_dir/$sample_name/${sample_name}_AllReads.sam

            samtools view -h $output_dir/$sample_name/${sample_name}_AllReads.sam > $output_dir/$sample_name/${sample_name}_AllReads.bam 
            grep "NH:i:1" $output_dir/$sample_name/${sample_name}_AllReads.sam > $output_dir/$sample_name/${sample_name}_uniqueReads.sam 
            samtools view -H $output_dir/$sample_name/${sample_name}_uniqueReads.sam  > $output_dir/$sample_name/${sample_name}_header.txt
            cat $output_dir/$sample_name/${sample_name}_header.txt $output_dir/$sample_name/${sample_name}_uniqueReads.sam > $output_dir/$sample_name/${sample_name}_uniqueReadsH.sam
            samtools view -Sb $output_dir/$sample_name/${sample_name}_uniqueReadsH.sam | samtools sort -o $output_dir/$sample_name/${sample_name}_unique_sorted.bam
            samtools index $output_dir/$sample_name/${sample_name}_unique_sorted.bam $output_dir/$sample_name/${sample_name}_unique_index.bam
            rm $output_dir/$sample_name/*.sam
            rm $output_dir/$sample_name/*_header.txt

            "$python" -m HTSeq.scripts.count \
            --format="$format" -a="$a" --order="$order" --stranded="$stranded" --mode="$mode" --nonunique="$nonunique" --type="$type" --idattr="$idattr" \
            $output_dir/$sample_name/${sample_name}_unique_index.bam "$GFF" > "$output_dir/$sample_name/Results.txt"

        done
    
    

    else
    
        local files=$(find "$fastq_dir" -type f \( -name "*_1.fastq" -o -name "*_1.fq" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \))

        for f in $files; do
            local sample_name=$(basename "$f" | sed -E 's/(\_1.fastq|\_1.fq|\_1.fastq\.gz|\_1.fq\.gz)$//')
            local dir_name=$(dirname "$f")
            
            if [[ $f == *".fq" ]]; then
            file1=$dir_name/$sample_name""_1.fq
            file2=$dir_name/${sample_name}_2.fq
            fi
            
            if [[ $f == *".fq.gz" ]]; then
            file1=$dir_name/$sample_name""_1.fq.gz
            file2=$dir_name/${sample_name}_2.fq.gz
            fi

            if [[ $f == *".fastq" ]]; then
            file1=$dir_name/$sample_name""_1.fastq
            file2=$dir_name/${sample_name}_2.fastq
            fi
            
            if [[ $f == *".fastq.gz" ]]; then
            file1=$dir_name/$sample_name""_1.fastq.gz
            file2=$dir_name/${sample_name}_2.fastq.gz
            fi
            
            
            mkdir -p "$output_dir/$sample_name"

            hisat2 -p "$threads" \
            -x $index_dir \
            --summary-file $output_dir/$sample_name/${f}_summary.txt \
            -1 $file1 \
            -2 $file2 \
            -S $output_dir/$sample_name/${sample_name}_AllReads.sam

            samtools view -h $output_dir/$sample_name/${sample_name}_AllReads.sam > $output_dir/$sample_name/${sample_name}_AllReads.bam
            grep "NH:i:1" $output_dir/$sample_name/${sample_name}_AllReads.sam > $output_dir/$sample_name/${sample_name}_uniqueReads.sam 
            samtools view -H $output_dir/$sample_name/${sample_name}_uniqueReads.sam  > $output_dir/$sample_name/${sample_name}_header.txt
            cat $output_dir/$sample_name/${sample_name}_header.txt $output_dir/$sample_name/${sample_name}_uniqueReads.sam > $output_dir/$sample_name/${sample_name}_uniqueReadsH.sam
            samtools view -Sb $output_dir/$sample_name/${sample_name}_uniqueReadsH.sam | samtools sort -o $output_dir/$sample_name/${sample_name}_unique_sorted.bam
            samtools index $output_dir/$sample_name/${sample_name}_unique_sorted.bam $output_dir/$sample_name/${sample_name}_unique_index.bam
            rm $output_dir/$sample_name/*.sam
            rm $output_dir/$sample_name/*_header.txt


             "$python" -m HTSeq.scripts.count \
            --format="$format" -a="$a" --order="$order" --stranded="$stranded" --mode="$mode" --nonunique="$nonunique" --type="$type" --idattr="$idattr" \
            $output_dir/$sample_name/${sample_name}_unique_index.bam "$GFF" > "$output_dir/$sample_name/Results.txt"
        done 
    fi 

echo "Complete"

}




