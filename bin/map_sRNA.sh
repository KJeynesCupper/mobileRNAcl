#!/bin/bash

# Function to perform alignment and clustering
map_sRNA() {
    local fasta_dir=""
    local index_dir=""
    local fastq_dir=""
    local output_dir=""
    local threads=6  
    local pad=200    
    local mincov=0.5 
    local dicermin=20
    local dicermax=24
    local dn_mirna="False"

    # Part 1 
   
    # Parse command-line options
    while getopts "FASTA:index:i:o:threads:pad:mincov:dicermin:dicermax:dn_mirna:" opt; do
        case ${opt} in
            FASTA) fasta_dir="$OPTARG" ;;   
            index) index_dir="$OPTARG" ;;    
            i) fastq_dir="$OPTARG" ;;    
            o) output_dir="$OPTARG" ;;   
            threads) threads="$OPTARG" ;;      
            pad) pad="$OPTARG" ;;          
            mincov) mincov="$OPTARG" ;;
            dicermin) dicermin="$OPTARG" ;;
            dicermax) dicermax="$OPTARG" ;; 
            dn_mirna) dn_mirna="$OPTARG" ;;      
            *) echo "Invalid option"; exit 1 ;;
        esac
    done

    # Check if all required arguments are provided
    if [[ -z "$fasta_dir" || -z "$index_dir" || -z "$fastq_dir" || -z "$output_dir" ]]; then
        echo "Usage: map_sRNA -r /path/to/fasta/dir -i /path/to/index/dir -f /path/to/fastq/dir -o /path/to/output/dir [-t threads] [-p pad] [-m mincov] [-a dicermin] [-b dicermax]"
        exit 1
    fi

    # Check if Bowtie index exists, and build if necessary
    if [[ ! -f "${index_dir}.1.ebwt" ]]; then
        echo "Bowtie index not found, building index..."
        bowtie-build --threads "$threads" "$fasta_dir" "$index_dir"
    fi

    # Create necessary directories
    mkdir -p "$output_dir/1_de_novo_detection/1_alignment"

    echo "Begining alignment of sequencing reads..."
    # List all FASTQ files in the directory
    local files=$(find "$fastq_dir" -type f \( -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" \))

    # Create or clear stats file
    echo -n > "$output_dir/1_de_novo_detection/stats_alignment.txt"

    # Alignment and clustering loop
    for f in $files; do
        # Extract the sample name by removing extensions
        local sample_name=$(basename "$f" | sed -E 's/(\.fastq|\.fq|\.fastq\.gz|\.fq\.gz)$//')

        echo "Mapping with Bowtie ... $sample_name" >> "$output_dir/1_de_novo_detection/stats_alignment.txt"

        # Perform mapping and sorting
        (bowtie -v 0 -a -m 1 --sam -x "$index_dir" "$f" | \
        samtools view -bS | \
        samtools sort -o "$output_dir/1_de_novo_detection/1_alignment/${sample_name}.bam") 2>> "$output_dir/1_de_novo_detection/stats_alignment.txt"

        echo "Clustering with ShortStack ... $sample_name" >> "$output_dir/1_de_novo_detection/stats_alignment.txt"

        # Cluster analysis with ShortStack
        ShortStack \
        --bamfile "$output_dir/1_de_novo_detection/1_alignment/${sample_name}.bam" \
        --genomefile "$fasta_dir" \
        --threads "$threads" \
        --pad "$pad" \
        --mincov "$mincov" \
        --nohp \
        --dicermin "$dicermin" \
        --dicermax "$dicermax" \
        --outdir "$output_dir/1_de_novo_detection/$sample_name" >> "$output_dir/1_de_novo_detection/stats_alignment.txt" 2>&1
    done

    # part 2 
    echo "Generating loci file..."

    # Define the merged output file path
    local locifile="$output_dir/locifile.txt"

    # Find all the Results.gff3 files and concatenate them
    find "$output_dir/1_de_novo_detection" -type f -name "Results.gff3" -exec cat {} + > "$locifile"

    # Remove duplicate lines
    sort "$locifile" | uniq > "${locifile}.tmp" && mv "${locifile}.tmp" "$locifile"

    echo "Loci file generate, and saved to: $locifile"
    

    # part 3 
    echo "Begining final clustering analysis..."
     # Create necessary directories
    mkdir -p "$output_dir/2_sRNA_results"

    echo "Begining alignment of sequencing reads..."

    local alignmentfiles="$output_dir/1_de_novo_detection/1_alignment/*.bam"

    if [[ "$dn_mirna" == "TRUE" ]]; then
    
    # clustering loop
    for f in $alignmentfiles; do
    f=${f##*/}
    sample_name=${f%.bam}

    
        # Cluster analysis with ShortStack
        ShortStack \
        --bamfile "$output_dir/1_de_novo_detection/1_alignment/${sample_name}.bam" \
        --genomefile "$fasta_dir" \
        --threads "$threads" \
        --pad "$pad" \
        --mincov "$mincov" \
        --nohp \
        --dn_mirna \
        --dicermin "$dicermin" \
        --dicermax "$dicermax" \
        --outdir "$output_dir/2_sRNA_results/$sample_name" >> "$output_dir/2_sRNA_results/stats_alignment.txt" 2>&1
    done

    else
    # clustering loop
    for f in $alignmentfiles; do
    f=${f##*/}
    sample_name=${f%.bam}

    
        # Cluster analysis with ShortStack
        ShortStack \
        --bamfile "$output_dir/1_de_novo_detection/1_alignment/${sample_name}.bam" \
        --genomefile "$fasta_dir" \
        --threads "$threads" \
        --pad "$pad" \
        --mincov "$mincov" \
        --nohp \
        --dicermin "$dicermin" \
        --dicermax "$dicermax" \
        --outdir "$output_dir/2_sRNA_results/$sample_name" >> "$output_dir/2_sRNA_results/stats_alignment.txt" 2>&1
    done
    fi

    

echo "Complete."
}

