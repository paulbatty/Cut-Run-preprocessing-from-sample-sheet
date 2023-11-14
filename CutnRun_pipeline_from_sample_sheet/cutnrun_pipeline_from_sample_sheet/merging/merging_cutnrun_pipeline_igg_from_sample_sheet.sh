#!/bin/bash

# Define the base directory. This is where the pipeline output will go
base_directory=/inpath/

# Define Slurm job settings
#SBATCH --output="$base_directory/CnR_pipe_%j.log"
#SBATCH --job-name=CnR
#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=8000
#SBATCH --error="$base_directory/CnR_pipe_%j.err"
#SBATCH --mail-type=end
#SBATCH --mail-user=email

# Exit script if any command fails
set -e

# Load necessary modules
module load Trim_Galore/0.6.6-GCC-10.2.0-Python-2.7.18
module load Bowtie2/2.4.2-GCC-10.2.0
module load SAMtools/1.15.1-GCC-11.2.0
module load BEDTools/2.30.0-GCC-10.3.0
module load deepTools/3.5.1-foss-2021a-Python-3.9.5

# Specify the input directory for BAM files
input_bam_directory="$base_directory/out/bam"

# Specify the output directory
output_bam_merged_directory="$base_directory/out/merged_bam"
output_bw_merged_directory="$base_directory/out/merged_bw"

# Specify the directory for logs and errors
merged_log_directory="$base_directory/merged_logs"
merged_error_directory="$base_directory/merged_errors"
fastqc_directory="$base_directory/fastqc"
trimming_directory="$base_directory/trimming"

# Create output directories if they do not exist
mkdir -p "$output_bam_merged_directory"
mkdir -p "$output_bw_merged_directory"
mkdir -p "$merged_log_directory"
mkdir -p "$merged_error_directory"
mkdir -p "$fastqc_directory"
mkdir -p "$trimming_directory"

# Read replicate pairs from the CSV file
replicate_input_path="$base_directory/replicate_pairs.csv"
IFS=','

# Debugging: Print the replicate pairs
echo "Replicate Pairs: $replicate_pairs"

# Get a timestamp for the log file
timestamp=$(date "+%Y%m%d_%H%M%S")

# Loop through each replicate pair
while read -r replicate1 replicate2
do
    # Remove newlines or other whitespace characters from the replicates
    replicate1_cleaned=$(echo "$replicate1" | tr -d '[:space:]')
    replicate2_cleaned=$(echo "$replicate2" | tr -d '[:space:]')

    # Use the input BAM directory to construct the paths to the input BAM files for each replicate
    input_bam_replicate1="${input_bam_directory}/${replicate1_cleaned}.bam"
    input_bam_replicate2="${input_bam_directory}/${replicate2_cleaned}.bam"

    # Construct the merged BAM file name by using the basename of replicate1 without the "-R1" suffix
    merged_basename=$(echo "$replicate1_cleaned" | sed 's/-R1//')
    merged_bam="${output_bam_merged_directory}/${merged_basename}_merged.bam"

    # Construct the full paths to the output files for the merged BAM file
    output_bam_rm_bklist="${output_bam_merged_directory}/${merged_basename}_blacklist_removed.bam"
    output_bw="${output_bw_merged_directory}/${merged_basename}.bw"
    orig_bam_file="${input_bam_directory}/${merged_basename}.bam"
    namesorted_bam="${output_bam_merged_directory}/${merged_basename}_namesorted.bam"
    fixmate_bam="${output_bam_merged_directory}/${merged_basename}_namesorted_fixmate.bam"
    positionsorted_bam="${output_bam_merged_directory}/${merged_basename}_namesorted_fixmate_positionsorted.bam"
    dups_removed_bam="${output_bam_merged_directory}/${merged_basename}_namesorted_fixmate_positionsorted_dupsremoved.bam"

    # Submit a Slurm job for processing
    sbatch --partition=mediumq --qos=mediumq --nodes=1 --ntasks=1 --cpus-per-task=4 --time=2-00:00:00 --mem=8G --error="${merged_error_directory}/all_samples_${timestamp}.err" --output="${merged_log_directory}/all_samples_${timestamp}.log" --job-name=cnr5migg \
    --wrap="
    # Step 1: Merge BAM files
    echo \"Processing merged barcode pair: $replicate1_cleaned-$replicate2_cleaned, Step 1: Merging BAM files started...\"
    samtools merge -f $merged_bam $input_bam_replicate1 $input_bam_replicate2 &&
    echo \"Processing merged barcode pair: $replicate1_cleaned-$replicate2_cleaned, Step 1: Merging BAM files finished.\"

    # Step 2: Sort BAM by read name
    echo \"Processing merged barcode pair: $replicate1_cleaned-$replicate2_cleaned, Step 2: Sorting BAM by read name...\"
    samtools sort -n $merged_bam -o $namesorted_bam &&
    echo \"Processing merged barcode pair: $replicate1_cleaned-$replicate2_cleaned, Step 2: Sorting BAM by read name finished.\"

    # Step 3: Fixmate information in BAM
    echo \"Processing merged barcode pair: $replicate1_cleaned-$replicate2_cleaned, Step 3: Fixing mate information...\"
    samtools fixmate -m $namesorted_bam $fixmate_bam &&
    echo \"Processing merged barcode pair: $replicate1_cleaned-$replicate2_cleaned, Step 3: Fixing mate information finished.\"

    # Step 4: Sort BAM by position
    echo \"Processing merged barcode pair: $replicate1_cleaned-$replicate2_cleaned, Step 4: Sorting BAM by position...\"
    samtools sort $fixmate_bam -o $positionsorted_bam &&
    echo \"Processing merged barcode pair: $replicate1_cleaned-$replicate2_cleaned, Step 4: Sorting BAM by position finished.\"

    # Step 5: Mark and remove duplicates from IgG control
    echo \"Processing merged barcode pair: $replicate1_cleaned-$replicate2_cleaned, Step 5: Marking and removing duplicates...\"
    samtools markdup -r $positionsorted_bam $dups_removed_bam &&
    echo \"Processing merged barcode pair: $replicate1_cleaned-$replicate2_cleaned, Step 5: Marking and removing duplicates finished.\"

    # Step 6: Remove reads overlapping with blacklist regions
    echo \"Processing merged barcode pair: $replicate1_cleaned-$replicate2_cleaned, Step 6: Removing reads overlapping with blacklist regions...\"
    intersectBed -abam $dups_removed_bam -b /Blacklist/CUTRun_Blacklist_hg38.bed -v > $output_bam_rm_bklist &&
    echo \"Processing merged barcode pair: $replicate1_cleaned-$replicate2_cleaned, Step 6: Removing reads overlapping with blacklist regions finished.\"

    # Step 7: Index the final BAM
    echo \"Processing merged barcode pair: $replicate1_cleaned-$replicate2_cleaned, Step 7: Indexing the final BAM...\"
    samtools index $output_bam_rm_bklist &&
    echo \"Processing merged barcode pair: $replicate1_cleaned-$replicate2_cleaned, Step 7: Indexing the final BAM finished.\"

    # Step 8: Generate coverage tracks
    echo \"Processing merged barcode pair: $replicate1_cleaned-$replicate2_cleaned, Step 8: Generating coverage tracks...\"
    bamCoverage -b $output_bam_rm_bklist -o $output_bw --normalizeUsing CPM --effectiveGenomeSize 2913022398 -p 4 --smoothLength 60 --binSize 10 &&
    echo \"Processing merged barcode pair: $replicate1_cleaned-$replicate2_cleaned, Step 8: Generating coverage tracks finished.\"

    # Clean up intermediate files
    echo \"Processing merged barcode pair: $replicate1_cleaned-$replicate2_cleaned, Cleaning up intermediate files...\"
    rm $namesorted_bam
    rm $fixmate_bam
    rm $positionsorted_bam
    echo \"Processing merged barcode pair: $replicate1_cleaned-$replicate2_cleaned, Clean up finished.\"
    echo \"All steps completed for merged barcode pair: $replicate1_cleaned-$replicate2_cleaned.\""
done < $replicate_input_path
