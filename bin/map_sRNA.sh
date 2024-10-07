#!/bin/bash

# Function to perform alignment and clustering
function map_sRNA() {

    local threads=6  
    local pad=200    
    local mincov=0.5 
    local dicermin=20
    local dicermax=24
    local dn_mirna=FALSE


while getopts ":f:i:o:x:t:p:m:d:c:n:" opt; do
            case ${opt} in
                f) FASTA="$OPTARG" ;;
                x) index="$OPTARG" ;;
                i) i="$OPTARG" ;;
                o) o="$OPTARG" ;;
                t) threads="$OPTARG" ;;
                p) pad="$OPTARG" ;;
                m) mincov="$OPTARG" ;;
                d) dicermin="$OPTARG" ;;
                c) dicermax="$OPTARG" ;;
                n) dn_mirna="$OPTARG" ;;
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
echo "FASTA: $FASTA"
echo "Input: $i"
echo "Output: $o"
echo "Index: $index"
echo "Threads: $threads"
echo "Pad: $pad"
echo "MinCov: $mincov"
echo "DicerMin: $dicermin"
echo "DicerMax: $dicermax"
echo "Dicer miRNA: $dn_mirna"

        
    # Check if Bowtie index exists, and build if necessary
    if [[ ! -f "${index}.1.ebwt" ]]; then
        echo "Bowtie index not found, building index..."
        bowtie-build --threads $threads "$FASTA" "$index"
    fi

    # Create necessary directories
    mkdir -p "$o/1_de_novo_detection/1_alignment"

    echo "Beginning alignment of sequencing reads..."
    # List all FASTQ files in the directory
    local files=$(find "$i" -type f \( -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" \))

    # Create or clear stats file
    echo -n > "$o/1_de_novo_detection/stats_alignment.txt"

    # Alignment and clustering loop
    for f in $files; do
        # Extract the sample name by removing extensions
        local sample_name=$(basename "$f" | sed -E 's/(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$//')

        echo "Mapping with Bowtie ... $sample_name" >> "$o/1_de_novo_detection/stats_alignment.txt"
        # Perform mapping and sorting
        (bowtie -v 0 -a -m 1 --sam "$index" "$f" | \
        samtools view -bS | \
        samtools sort -o "$o/1_de_novo_detection/1_alignment/${sample_name}.bam") 2>> "$o/1_de_novo_detection/stats_alignment.txt"

        echo "Clustering with ShortStack ... $sample_name" >> "$o/1_de_novo_detection/stats_alignment.txt"

        # Cluster analysis with ShortStack
        ShortStack \
        --bamfile "$o/1_de_novo_detection/1_alignment/${sample_name}.bam" \
        --genomefile "$FASTA" \
        --threads "$threads" \
        --pad "$pad" \
        --mincov "$mincov" \
        --nohp \
        --dicermin "$dicermin" \
        --dicermax "$dicermax" \
        --outdir "$o/1_de_novo_detection/$sample_name" >> "$o/1_de_novo_detection/stats_alignment.txt" 2>&1
    done

    # Part 2 
    echo "Generating loci file..."

    # Define the merged output file path
    local locifile="$o/locifile.txt"

    # Find all the Results.gff3 files and concatenate them
    find "$o/1_de_novo_detection" -type f -name "Results.gff3" -exec cat {} + > "$locifile"

    # Remove duplicate lines
    sort "$locifile" | uniq > "${locifile}.tmp" && mv "${locifile}.tmp" "$locifile"

    echo "Loci file generated, and saved to: $locifile"
    

    # Part 3 
    echo "Beginning final clustering analysis..."
    # Create necessary directories
    mkdir -p "$o/2_sRNA_results"


    local alignmentfiles="$o/1_de_novo_detection/1_alignment/*.bam"

    if [[ "$dn_mirna" == "TRUE" ]]; then
    
    # Clustering loop
    for f in $alignmentfiles; do
        f=${f##*/}
        sample_name=${f%.bam}

        # Cluster analysis with ShortStack
        ShortStack \
        --bamfile "$o/1_de_novo_detection/1_alignment/${sample_name}.bam" \
        --genomefile "$FASTA" \
        --threads "$threads" \
        --pad "$pad" \
        --mincov "$mincov" \
        --nohp \
        --dn_mirna \
        --dicermin "$dicermin" \
        --dicermax "$dicermax" \
        --outdir "$o/2_sRNA_results/$sample_name" >> "$o/2_sRNA_results/stats_alignment.txt" 2>&1
    done

    else
        # Clustering loop
        for f in $alignmentfiles; do
            f=${f##*/}
            sample_name=${f%.bam}

            # Cluster analysis with ShortStack
            ShortStack \
            --bamfile "$o/1_de_novo_detection/1_alignment/${sample_name}.bam" \
            --genomefile "$FASTA" \
            --threads "$threads" \
            --pad "$pad" \
            --mincov "$mincov" \
            --nohp \
            --dicermin "$dicermin" \
            --dicermax "$dicermax" \
            --outdir "$o/2_sRNA_results/$sample_name" >> "$o/2_sRNA_results/stats_alignment.txt" 2>&1
        done
    fi

    echo "Complete."

}
