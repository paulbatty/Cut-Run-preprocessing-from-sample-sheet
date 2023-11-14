#!/bin/bash

# Define the base directory
base_directory="/inpath/"

# Define Slurm job settings
#SBATCH --output="$base_directory/bamcov_bdg_SEACR/bamcov_bdg_%j.log"
#SBATCH --job-name=bamcovbdg
#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=8000
#SBATCH --error="$base_directory/bamcov_bdg_SEACR/bamcov_bdg_%j.err"
#SBATCH --mail-type=end
#SBATCH --mail-user=email

# Load necessary modules
module load deepTools/3.5.1-foss-2021a-Python-3.9.5

# Specify the input directory for BAM files
input_bam_directory="$base_directory/out/merged_bam"

# Specify the output directory for bedGraph files
output_bedgraph_directory="$base_directory/out/merged_bedgraph"

# Specify the directory for logs and errors
bedgraph_log_directory="$base_directory/log/merged_bedgraph_logs"
bedgraph_error_directory="$base_directory/errors/merged_bedgraph_errors"

# Create output directories if they do not exist
mkdir -p "$output_bedgraph_directory"
mkdir -p "$bedgraph_log_directory"
mkdir -p "$bedgraph_error_directory"

# Get a timestamp for the log and error file
timestamp=$(date "+%Y%m%d_%H%M%S")

# Define the log and error files with timestamp
log_file="$bedgraph_log_directory/bamcov_bdg_merged_${timestamp}.log"
error_file="$bedgraph_error_directory/bamcov_bdg_merged_${timestamp}.err"

# Loop through each BAM file ending with *blacklist_removed.bam
for i in "$input_bam_directory"/*blacklist_removed.bam
do
    NAME=`basename $i | cut -d"_" -f1`
    sbatch --partition=mediumq --qos=mediumq --nodes=1 --ntasks=1 --cpus-per-task=1 --time=2-00:00:00 --mem=8G --error=$error_file --output=$log_file --job-name=bdgcnr4m --wrap="bamCoverage -b $i -o $output_bedgraph_directory/${NAME}.bdg -of bedgraph --normalizeUsing CPM --effectiveGenomeSize 2913022398 -p 12 --smoothLength 60  --binSize 10"
done
