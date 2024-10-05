#!/bin/bash
#SBATCH --qos bbdefault
#SBATCH --ntasks 40
#SBATCH --nodes 1
#SBATCH --time 100

module purge; module load bluebear
module load SAMtools/1.10-GCC-8.3.0

# ################################################################################
# ### 1 Raw reads
# ################################################################################
# for file in ./1_raw/*.fq.gz; do
#     read_count=$(wc -l < "$file")
#     # Calculate the number of reads by dividing the total lines by 4 (each read has 4 lines)
#     read_count=$((read_count / 4))
#     echo "$file $read_count" >> ./6_mapping_metrics/1_raw_read_metrics.txt
# done
#
#
# ################################################################################
# ### 2 Trimmed reads
# ################################################################################
# for file in ./2_trimmed/*.fq.gz; do
#     read_count=$(wc -l < "$file")
#     # Calculate the number of reads by dividing the total lines by 4 (each read has 4 lines)
#     read_count=$((read_count / 4))
#     echo "$file $read_count" >> ./6_mapping_metrics/2_trimmed_read_metrics.txt
# done
#
#
#
# ################################################################################
# ### 3 mapped reads
# ################################################################################
#
# # Output file to store the results
# output_file="./6_mapping_metrics/3_mapping_metrics.txt"
#
# # Clear the output file if it already exists
# > "$output_file"
#
#
# # Loop through sample folders
# for sample_dir in ./3_de_novo_detection/1_alignment/*.bam; do
#         # Extract the sample name from the directory name
#         sample_name=$(basename "$sample_dir")
#         sample_name=${sample_name%.bam}
#
#         # BAM file path
#          bam_file=${sample_dir}
#
# ##### spikes
# # extract reads mapped to each genome
# samtools view -h  "$bam_file" | grep -E '^\@|A_' | samtools view -b -o "./6_mapping_metrics/$sample_name""_A.bam" -
# samtools view -h  "$bam_file" | grep -E '^\@|B_' | samtools view -b -o "./6_mapping_metrics/$sample_name""_B.bam" -
#
# # Number of mapped reads:
# total_mapped=$(samtools view -c -F 0x4 "$bam_file")
# unmapped=$(samtools view -c -f 4 "$bam_file")
# total_reads=$(samtools view -c "$bam_file")
#
# mapped_to_A_total=$(samtools view -c -F 0x4 "./6_mapping_metrics/$sample_name""_A.bam")
# mapped_to_B_total=$(samtools view -c -F 0x4 "./6_mapping_metrics/$sample_name""_B.bam")
#
#
# # Print the results to the output file
# echo "Sample: $sample_name" >> "$output_file"
# echo "Total reads: $total_reads" >>  "$output_file"
# echo "Mapped Reads: $total_mapped" >> "$output_file"
# echo "Unmapped Reads: $unmapped" >> "$output_file"
# echo "---" >> "$output_file"
# echo "Reads mapped to Genome A: $mapped_to_A_total" >> "$output_file"
# echo "Reads mapped to Genome B: $mapped_to_B_total" >> "$output_file"
# echo "---" >> "$output_file"
# echo "" >> "$output_file"  # Add a blank line for separation
#
# done
# rm ./6_mapping_metrics/*.bam

################################################################################
### As CVS file
################################################################################

# Output file to store the results
output_file="./6_mapping_metrics/3_mapping_metrics.csv"
echo "Sample, Raw_reads, Cleaned_Reads, Total_reads,Mapped_Reads, Unmapped_Reads, Mapped_to_GenomeA, Mapped_to_GenomeB" > "$output_file"; \

for file in ./1_raw/*.fq.gz; do

  sample_name=$(basename "$file")
  sample_name=${sample_name%.fq.gz}

# 1 - caclculate raw reads
rawread_count=$(wc -l < "$file")
rawread_count=$((rawread_count / 4))

# 2 - trim count
trim_file="./2_trimmed/$sample_name""_trimmed.fq.gz"
trimread_count=$(wc -l < "$trim_file")
trimread_count=$((trimread_count / 4))


# 3 - mapping metrics
alignment_file="./3_de_novo_detection/1_alignment/$sample_name"".bam"

# extract reads mapped to each genome
samtools view -h  "$alignment_file" | grep -E '^\@|A_' | samtools view -b -o "./6_mapping_metrics/$sample_name""_A.bam" -
samtools view -h  "$alignment_file" | grep -E '^\@|B_' | samtools view -b -o "./6_mapping_metrics/$sample_name""_B.bam" -

# Number of mapped reads:
total_mapped=$(samtools view -c -F 0x4 "$alignment_file")
unmapped=$(samtools view -c -f 4 "$alignment_file")
total_reads=$(samtools view -c "$alignment_file")

mapped_to_A_total=$(samtools view -c -F 0x4 "./6_mapping_metrics/$sample_name""_A.bam")
mapped_to_B_total=$(samtools view -c -F 0x4 "./6_mapping_metrics/$sample_name""_B.bam")


# Print the results to the output file
echo "$sample_name, $rawread_count, $trimread_count, $total_reads,$total_mapped, $unmapped, $mapped_to_A_total, $mapped_to_B_total "  >> "$output_file"
rm ./6_mapping_metrics/*.bam


done
