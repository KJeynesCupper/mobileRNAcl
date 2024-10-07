#!/bin/bash

function RNAmergeGenomes() {
    # Initialize an array for input files
    input_files=()
    
    # Parse the input options
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -i) # When -i is encountered, shift to the next argument and gather all subsequent arguments as input files until -o is encountered
                shift
                while [[ "$1" != "-o" && $# -gt 0 ]]; do
                    input_files+=("$1")
                    shift
                done
                ;;
            -o) # Output file
                output_file="$2"
                shift 2
                ;;
            *) # Catch any unrecognized option or argument
                echo "Usage: $0 -i <input_files>... -o <output_file>"
                return 1
                ;;
        esac
    done

    # Check point: are there both input_files and output_file 
    if [[ ${#input_files[@]} -eq 0 || -z "$output_file" ]]; then
        echo "You must provide at least one input file with -i and an output file with -o."
        return 1
    fi

    # Temporary file to store the concatenated result
    temp_file=$(mktemp)

    # Iterate over the input files
    for i in "${!input_files[@]}"; do
        file="${input_files[$i]}"
        prefix=$(printf "%b_" $(printf '\\x%X' $((65 + i)))) #  add prefixes to pseudomolecule. 
        
        if [[ -f "$file" ]]; then
            if [[ "$file" == *.gz ]]; then
                # for  gzipped files
                gunzip -dc "$file" | sed "s/^>/>$prefix/" >> "$temp_file"
            
            else
                # Process regular files
                sed "s/^>/>$prefix/" "$file" >> "$temp_file"
            fi
        else
            echo "File not found: $file"
            return 1
        fi
    done

    
    # if user wants tozip output 
    if [[ "$output_file" == *.gz ]]; then
        gzip -c "$temp_file" > "$output_file"
    else
        mv "$temp_file" "$output_file"
    fi
    echo "Genome reference files are merged. Merged Reference file is $output_file"
}
