# integrato
Mapping the integration sites from UDiTaS experiments with the eePASSIGE system

Overview

Integrato is a tool for identifying and quantifying the on- and off-target integration events of programmable gene integration technologies. Data are collected with a UDiTaS experiment, then analyzed with the script here.  
Summary of the UDiTaS experiment. Treatment of cells with integrase results in on-target integration of the donor DNA as well as off-target integration at unknown loci. Tagmentation with Tn5* introduces Illumina adapters randomly across the genome, providing handles for PCR amplification and sequencing of on- and off-target sites. 

![image](https://github.com/Nicholas-Krasnow/integrato/assets/119907096/cd5fb2ff-41df-48e3-ad8a-c6734bf7f7d5)
 
Summary of UDiTaS data analysis pipeline. Sequencing reads generated in the experiment are first processed with quality filtering as well as filtering for reads with the donor-specific primer and att subsequence. Reads are then deduplicated by the UMI sequence in the donor DNA. Filtered and deduplicated reads are aligned to the genome with BWA-MEM and clustered into neighborhoods with Bedtools which allows for the identification and quantification of genomic integration sites. 

Experiment

The UDiTaS experiment and data collection should be conducted exactly as described in the paper with the integrase of your choosing. When naming samples for the Miseq run, start the sample name with the target site followed by a dash, and don’t use any underscores.
e.g., AAVS1-WT-Bxb1-2weeks

Analyzing the sequencing data

The resulting fastqs from the experiment can be analyzed with the following steps on MacOS. 

Install system requirements:
1. Anaconda: https://www.anaconda.com/download
2. Homebrew: https://brew.sh/
3. bwa: run the command in the terminal ```brew install bwa```  
4. bedtools: run the command in the terminal brew install bedtools 

Setup environment: This step will create a conda virtual environment with the necessary Python installation and necessary packages. Only need to do this once.

1.	Navigate to the location in your file system where you want to run the analysis (e.g., ~/Documents)
2.	Download the entire repo to this location
3.	Open a terminal window and change directory to this one with “cd” command (e.g., cd ~/Documents/integrato)
4.	Run the environment setup script with the following terminal command:

sh create_env.sh

enter “y” when prompted.

Run the analysis: the data can now be analyzed with the few remaining steps. Start here if the conda environment has already been set up.

1.	In the terminal, change directory to the integrato analysis folder where the repo code was downloaded (e.g., cd ~/Documents/integrato)
2.	Copy the miseq output R1 fastqs from your experiment into this folder. Files should be compressed (.gz) and do not include the index files.
3.	Activate the conda environment with the terminal command:

conda activate integrato_env

4.	Run the analysis code with the following terminal command:

sh runISM_0.sh 

This command will execute all the code to analyze the reads and generate output quantification files.

Note: the code is written for identifying integrations of the donor plasmid and donor-specific primer used in the paper. Analysis with donor plasmids will require changes to the source code and a different indexed plasmidref reference. 


Interpreting results
Running the code should have generated a folder for each target site (e.g., ‘AAVS1’). In this folder, the subfolder “output_clustered_multiple_reads”, which contains the output files quantifying all clustered integration sites that were mapped with at least two reads. When the aligned genomic sequence is found directly adjacent to the att subsequence from the donor plasmid (allowing up to one mismatch from sequencing error), the expected pre-integration att sequence on the genome is provided. Note: UDiTaS experiments can produce false-positive integration loci, such as alignments arising from integration into the pegRNA plasmids. Nominated sites should be manually inspected, further filtered for known artifacts, and validated with an orthogonal assay (ddPCR).
