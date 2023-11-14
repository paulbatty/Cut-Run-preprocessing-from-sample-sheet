#!/bin/bash

# Define the base directory. This is where the pipeline output will go
base_directory=/inpath/

# Define Slurm job settings
#SBATCH --output=$base_directory/CnR_pipe_%j.log
#SBATCH --job-name=CnR
#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=8000
#SBATCH --error $base_directory/CnR_pipe_%j.err
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

# Specify the input directory. This is where the FASTQ files are
input_directory=/inpath/files/

# Specify the output directory
output_directory=$base_directory/out/bigwig
output_bam_directory=$base_directory/out/bam

# Specify the output directory for bedGraph files
output_bedgraph_directory=$base_directory/out/bedgraph

# Specify the directory for logs and errors
log_directory=$base_directory/log
error_directory=$base_directory/errors
fastqc_directory=$base_directory/fastqc
trimming_directory=$base_directory/trimming

# Create output directories if they do not exist
mkdir -p "$output_directory"
mkdir -p "$output_bam_directory"
mkdir -p "$output_bedgraph_directory"
mkdir -p "$log_directory"
mkdir -p "$error_directory"
mkdir -p "$fastqc_directory"
mkdir -p "$trimming_directory"

# Read barcodes from sample sheet
barcode_input_path=$base_directory/sample_sheet.csv
barcodes=$(cat "$barcode_input_path")

# Debugging: Print the barcodes
echo "Barcodes: $barcodes"

# Get a timestamp for the log file
timestamp=$(date "+%Y%m%d_%H%M%S")

# Loop through each barcode
for barcode in $barcodes
do
    # Remove newlines or other whitespace characters from the barcode
    barcode_cleaned=$(echo "$barcode" | tr -d '[:space:]')
    
    # Get the input file names for the barcode
    file_r1="${input_directory}/${barcode_cleaned}_R1_001.fastq.gz"
    file_r2="${input_directory}/${barcode_cleaned}_R2_001.fastq.gz"
    
    # Get the trimmed file names
    trimmed_r1="${trimming_directory}/${barcode_cleaned}_R1_001_val_1.fq.gz"
    trimmed_r2="${trimming_directory}/${barcode_cleaned}_R2_001_val_2.fq.gz"

    # Construct the full paths to the output files
    # Define the final BAM file paths
    output_bam_rm_bklist="${output_bam_directory}/${barcode_cleaned}_blacklist_removed.bam"
    output_bw="${output_directory}/${barcode_cleaned}.bw"
    output_bedgraph="${output_bedgraph_directory}/${barcode_cleaned}.bdg"
    sam_file="${trimming_directory}/${barcode_cleaned}.sam"
    bam_file="${output_bam_directory}/${barcode_cleaned}.bam"
    namesorted_bam="${output_bam_directory}/${barcode_cleaned}_namesorted.bam"
    fixmate_bam="${output_bam_directory}/${barcode_cleaned}_namesorted_fixmate.bam"
    positionsorted_bam="${output_bam_directory}/${barcode_cleaned}_namesorted_fixmate_positionsorted.bam"
    dups_removed_bam="${output_bam_directory}/${barcode_cleaned}_namesorted_fixmate_positionsorted_PCR_dup_rm.bam"

    # Submit a Slurm job for processing
    sbatch --partition=mediumq --qos=mediumq --nodes=1 --ntasks=1 --cpus-per-task=4 --time=2-00:00:00 --mem=8G --error="${error_directory}/all_samples_${timestamp}.err" --output="${log_directory}/all_samples_${timestamp}.log" --job-name=cnr4igg \
    --wrap="
    # Step 1: Trim and preprocess data
    echo \"Processing barcode: $barcode_cleaned, Step 1: Trimming and preprocessing started...\"
    trim_galore --illumina --fastqc --paired $file_r1 $file_r2 -o $trimming_directory &&
    echo \"Processing barcode: $barcode_cleaned, Step 1: Trimming and preprocessing finished.\"

    ## Step 2: Align reads using Bowtie2
    echo \"Processing barcode: $barcode_cleaned, Step 2: Alignment started...\"
    bowtie2 -p 4 --local --very-sensitive --no-mixed --no-discordant --dovetail -1 $trimmed_r1 -2 $trimmed_r2 -x /hg38/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome -S $sam_file &&
    echo \"Processing barcode: $barcode_cleaned, Step 2: Alignment finished.\"

    # Step 3: Convert SAM to BAM
    echo \"Processing barcode: $barcode_cleaned, Step 3: Converting SAM to BAM...\"
    samtools view -bS $sam_file -o $bam_file &&
    echo \"Processing barcode: $barcode_cleaned, Step 3: Converting SAM to BAM finished.\"

    # Step 4: Sort BAM by read name
    echo \"Processing barcode: $barcode_cleaned, Step 4: Sorting BAM by read name...\"
    samtools sort -n $bam_file -o $namesorted_bam &&
    echo \"Processing barcode: $barcode_cleaned, Step 4: Sorting BAM by read name finished.\"

    # Step 5: Fixmate information in BAM
    echo \"Processing barcode: $barcode_cleaned, Step 5: Fixing mate information...\"
    samtools fixmate -m $namesorted_bam $fixmate_bam &&
    echo \"Processing barcode: $barcode_cleaned, Step 5: Fixing mate information finished.\"

    # Step 6: Sort BAM by position
    echo \"Processing barcode: $barcode_cleaned, Step 6: Sorting BAM by position...\"
    samtools sort $fixmate_bam -o $positionsorted_bam &&
    echo \"Processing barcode: $barcode_cleaned, Step 6: Sorting BAM by position finished.\"

    # Step 7: Mark duplicates in IgG control
    echo \"Processing barcode: $barcode_cleaned, Step 7: Marking and removing PCR duplicates...\"
    samtools markdup -r $positionsorted_bam $dups_removed_bam &&
    echo \"Processing barcode: $barcode_cleaned, Step 7: Marking and removing PCR duplicates finished.\"

    # Step 8: Remove reads overlapping with blacklist regions
    echo \"Processing barcode: $barcode_cleaned, Step 8: Removing reads overlapping with blacklist regions...\"
    intersectBed -abam $dups_removed_bam -b /hg38/Sequence/Blacklist/CUTRun_Blacklist_hg38.bed -v > ${output_bam_rm_bklist} &&
    echo \"Processing barcode: $barcode_cleaned, Step 8: Removing reads overlapping with blacklist regions finished.\"

    # Step 9: Index the final BAM
    echo \"Processing barcode: $barcode_cleaned, Step 9: Indexing the final BAM...\"
    samtools index ${output_bam_rm_bklist} &&
    echo \"Processing barcode: $barcode_cleaned, Step 9: Indexing the final BAM finished.\"

    # Step 10: Generate coverage tracks
    echo \"Processing barcode: $barcode_cleaned, Step 10: Generating coverage tracks...\"
    bamCoverage -b ${output_bam_rm_bklist} -o $output_bw --normalizeUsing CPM --effectiveGenomeSize 2913022398 -p 4 --smoothLength 60 --binSize 10 &&
    echo \"Processing barcode: $barcode_cleaned, Step 10: Generating coverage tracks finished.\"

    # Step 11: Generate bedGraph file
    echo \"Processing barcode: $barcode_cleaned, Step 11: Generating bedGraph file...\"
    bamCoverage -b $output_bam_rm_bklist -o $output_bedgraph --normalizeUsing CPM --effectiveGenomeSize 2913022398 -p 4 --smoothLength 60 --binSize 10 &&
    echo \"Processing barcode: $barcode_cleaned, Step 11: Generating bedGraph file finished.\"

    # Clean up intermediate files
    # Do not remove the bam file, you will use as the input file for merging
    echo \"Processing barcode: $barcode_cleaned, Cleaning up intermediate files...\"
    mv ${trimming_directory}/${barcode_cleaned}_R1_001_val_1.fq.gz_fastqc.html $fastqc_directory
    mv ${trimming_directory}/${barcode_cleaned}_R2_001_val_2.fq.gz_fastqc.html $fastqc_directory
    rm ${trimming_directory}/${barcode_cleaned}_R1_001_val_1.fq.gz_fastqc.zip
    rm ${trimming_directory}/${barcode_cleaned}_R2_001_val_2.fq.gz_fastqc.zip
    rm $sam_file
    rm $namesorted_bam
    rm $fixmate_bam
    rm $positionsorted_bam
    rm $dups_removed_bam
    rm $trimmed_r1
    rm $trimmed_r2
    echo \"Processing barcode: $barcode_cleaned, Clean up finished.\"
    echo \"All steps completed for barcode: $barcode_cleaned.\""
done
