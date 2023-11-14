# CutnRun_pipeline_from_sample_sheet

Pipeline written in bash that takes a sample sheet with parcodes and parses them to generate the output data. Relative paths are used.

# cutnrun_pipeline_from_sample_sheet >>> no merging
- The pipeline is designed to be run on a scientific cluster using slurm. A sample sheet is provided with the unique barcodes for each file before the file ending. The sample sheet is a csv file with a single column, 'barcodes', and the unique barcodes found in the input directory containing the FastQ files to be analysed.

- This sample sheet can be generated using the following iPython notebook: generate_sample_sheets.ipynb

- The notebook can be run locally following installation of anaconda and Visual Studio Code. The environment yml file used to run the notebook can also be found in the repository. Instructions for how to set up this environment can be found in the rd_me.text file.

- In brief the notebook extracts the file name from a directory of files. As there are two different pipelines for dealing with IgG and non-IgG samples, a regular expression is set up to identify all sample names that contain the word IgG. These are output as one CSV file. The other barcodes are output as a second CSV file. Additional regular expressions can be set up if there are for example multiple experiments in the same FastQ directory and you don't want to analyse everything at once.

- Each of the pipelines (what they do exactly is explained below) requires a base directory that you should specify in the nobackup directory on the CeMM cluster. Within this directory you should clone the latest version of the pipeline that you need and also copy the sample sheet for the samples you want to analyse into this folder.

Within the pipeline file, you should change the following parameters:
- Input path >>> This is the directory where the FastQ files are that you want to analyse.
- Base_directory >>> This is the directory where the pipeline is saved and where all the resultant directories with the output files will be generated.
- The job name within the batch command
- If appropriate, change the slurm requirements e.g. number of cores and cpus per tasks. This script has n=1 cores, n=1 tasks, and n=4 cpus per task.

- When you want to run the pipeline you should navigate to the directory when it is and run it as a bash command in the working directory: e.g. >>> ./no_merge_cutnrun_pipeline_non_igg_sample_sheet.sh


All of the output directories will be automatically created. The following folders are created.
- A 'log' folder will be created with a timestamped log of all the samples. When a particular step has been finished for a barcode it will be printed to the log file.
- An 'errors' folder will do the same for any error messages.
- A 'trimming' folder is output where the aligned fastqc.gz files after trimming are stored. These are the input for generating the sam file.
- A 'fastqc' folder where all the fastqc.html files will be stored.
- An 'out' folder with 'bam' and 'bigwig' subfolders where the bam and bigwig files will be stored.
- The barcodes will be passed to the script and iterated over in order to perform the different steps of the pipeline. A single log file is generated that comments on when a step of the pipeline has been finished for a given barcode.

# cutnrun_pipeline_from_sample_sheet >>> no merging, step by step
There are two main pipelines, which take the csv file with IgG or non-IgG barcodes as an input. In the pipeline without merging, only single replicates are processed. The difference between the IgG and non-IgG barcodes pipeline is that for IgG samples, there is an additional processing step to remove PCR duplicates before removing blacklisted regions. This step is absent for the sample pipeline.

- Step 1: Trim and preprocess data using trim_galore

- Step 2: Align reads to the reference genome using Bowtie2

- Step 3: Convert SAM to BAM

- Step 4: Sort BAM by read name

- Step 5: Fixmate information in BAM

- Step 6: Sort BAM by genomic position

- ** Step 7 - IgG only ** removal of PCR duplicates

- Step 7: Remove reads overlapping with blacklist regions, provided as a bed file and specific for CutNRun experiments

- Step 8: Index the final BAM (Creating BAM index)

- Step 9: Generate coverage tracks (Creating coverage tracks in bigWig format) from which you can call peaks

- Step 10: Removal of intermediate files that are of no use for the downstream analysis.

- By default the following files are removed, this can be changed as necessary:

>>> sam files, all intermediate bam files, fastqc.zip files, trimmed fastqc.gz files after alignment.

- The initial bam files are retained by default as these are the input for the pipeline for merging. The final processed bam files after removing blacklisted regions are also retained.

- The final bigwig files will be saved as an output in the pipeline directory.


# cutnrun_pipeline_from_sample_sheet >>> merging
- If you want to merge your individual replicates you should run this pipeline. It takes the bam files for each replicate as an input and runs through the rest of the pipeline as for the single replicates. **You should run this after first running the pipeline for single replicates as it requires the bam files from that as an input.**

- This also takes a sample sheet. You can run the ipynb script find_replicate_pairs_all.ipynb to parse the replicate for a given condition into a csv file titled 'replicate_pairs.csv' which has 'replicate1', and 'replicate2' as its headers. It runs a regular expression analysis to find file names which share a common substring in the fastq.gz files.  The ipynb script also prints the output so you can do a sanity check to make sure it is working correctly. Alternatively you can just check the csv file. **You must make sure that your replicates have unique common strings in order for this to work**

- The replicate_pairs csv file is then fed into the merging pipeline. The logic is the same as for the pipeline without merging. You need to specify a base directory, in which the bam files are stored and also specify the location of the csv file. The replicate_pairs.csv file should be put in the base directory. All output folders, log files, error files, and output files will be output automatically.

# cutnrun_pipeline_from_sample_sheet >>> Add bedgraphs
- The input for the peak calling using SEACR or MACS requires the coverage file in a bedgraph format rather than bigwig. Here you find pipelines that generate all the files necessary for the peak calling, for both single and merged replicates.

# generate_bedgraphs_from_bam_separate_script
- It may be that you want to generate the bedgraph files for peak calling independently of the main preprocessing pipeline. Here you find an example script that takes the final processed bam file as an input and outputs the bedgraph files. The bedgraph files can then be used to call peaks.
