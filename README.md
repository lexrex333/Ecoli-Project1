# Running this Pipeline

## Dependencies
Python modules:
* [csv](https://docs.python.org/3/library/csv.html)
* [subprocess](https://docs.python.org/3/library/subprocess.html)
* [os](https://docs.python.org/3/library/os.html)
* [pandas](https://pandas.pydata.org/)
* [itertools import combinations](https://docs.python.org/3/library/itertools.html)
* [numpy](https://numpy.org/)
* [matplotlib](https://matplotlib.org/stable/)
* [seaborn](https://github.com/mwaskom/seaborn)

Tools:
* [FastQC](https://github.com/s-andrews/FastQC)
* [MultiQC](https://github.com/MultiQC/MultiQC)
* [Trimmomatic](https://github.com/timflutre/trimmomatic)
* [SPAdes](https://github.com/ablab/spades)
* [FastANI](https://github.com/ParBLiSS/FastANI)
* [MASH](https://github.com/marbl/mash)

For downloading the sequencing reads:
* [Prefetch and Fasterq-Dump](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump)

## Data Retrieval 
* Download raw sample sequencing reads of 12 LTEE genomes from NCBI SRA (1 from each population)
* Retrieve ancestral strain genome [REL606](https://www.ncbi.nlm.nih.gov/nuccore/NC_012967.1)

If you wanted to download your own reads, you can use this script to do so:
* [Downloading Reads from NCBI](https://github.com/lexrex333/Ecoli-Project1/blob/main/Get_Fastqs.py)

## Input Data:
* 12 _E. coli_ genomes from the Long Term Evolution Experiment (LTEE)
* One ancestor (REL606) genome from the LTEE
## References 
* [Wiki Link](https://github.com/lexrex333/Ecoli-Project1/wiki) for project 
* [Introduction to the Long-Term Evolution Experiment (LTEE)](https://the-ltee.org/about/) 
* [Tempo and mode of genome evolution in a 50,000-generation experiment](https://www.nature.com/articles/nature18959)

## Installing Software Dependencies
Make sure all of the tools are installed within the directory where you will be running everything. 

In case you have issues with MultiQC, MASH, or Trimmomatic, please refer to this for additional help - also see their githubs linked up above. 

Installing MultiQC: 
```bash

git clone https://github.com/MultiQC/MultiQC.git

cd MultiQC

pip install .
```

Installing MASH:
```bash

wget https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar

tar -xvf mash-Linux64-v2.3.tar

cp mash-Linux64-v2.3/mash .
```

Installing Trimmomatic:
```bash
Micky add stuff here.
```
## What this Pipeline Does and What it Outputs
This pipeline is comparing 12 _E. coli_ strains to their ancestor strain, from which they diverged into 12 different populations. To compare these strains to one another, we used 3 different metrics that have not been used to compare these specific strains in other studies using the same strains. These metrics are ANI, MASH, and dDDH. 

Previously, we downloaded all of the SRR FASTQ files from NCBI and then ran them through Prefetch and Fasterq-Dump. This resulted in FASTQ files that could be used for running this pipeline. 

### 1. Discovering Read Quality and Trimming
The pipeline first starts off with using FastQC to access the quality of each of the raw reads. MultiQC then gives a single summarized file which compiles all of the output from FastQC. This makes it easier to read, rather than having a separate file for each sample. 

The reads are then trimmed using Trimmomatic. Then, these trimmed reads are run again through FastQC and MultiQC to see the differences in read quality compared to before they were trimmed. 

### 2. Assembling Genomes to Use for Comparisons
Using the trimmed reads, SPAdes creates assembled genomes that can then be used to compare strains to one another. This will output FASTA files that will be used in the following 3 metrics. 

### 3. Three Comparison Metrics
After obtaining the FASTA files, ANI, MASH, and dDDH will take each of the population members' FASTA files and compare them to the ancestor strain. This will give a total of 12 comparisons (each most recent population member compared to the ancestor strain). Each of these metrics will output a CSV or TSV file containing their results. 

### 4. Comparing Metrics
To compare the results from the 3 different metrics, a TSV file combining all 3 metrics' output information is created.

### 5. Visualization
Using the information from this TSV file, the pipeline creates a heatmap and a bar chart showing the different similarity calculations between each of the metrics.



# How to Run This Pipeline

### 1. Retrieving necessary files
In order for this pipeline to run, you must obtain the following files from this repository:

[Pipeline script](https://github.com/lexrex333/Ecoli-Project1/blob/main/pipeline.py)

[CSV of all genome information](https://github.com/lexrex333/Ecoli-Project1/blob/main/recent_strains.csv)

You may obtain these either by downloading them manually and uploading them to your workspace, or by cloning this repository and pulling down the files:
```
git clone https://github.com/lexrex333/Ecoli-Project1.git
```
Another item that is necessary to run this pipeline is a directory called "Fastqs" which contains all of the raw FASTQ read files to be analyzed.
This can be obtained in one of two ways:
#### If using full-length FASTQ reads:
Run [this script](https://github.com/lexrex333/Ecoli-Project1/blob/main/Get_Fastqs.py) to download genomes directly from NCBI.
#### If using shortened sample reads:
A folder containing the first 250,000 reads from each FASTQ can be found in this directory on the class server:
```
/home/project1/Ecoli-Project1/Fastqs
```
It is imperative that you copy this directory into your own workspace to avoid pathway-related errors when running the pipeline. You must leave the directory named "Fastqs".
```
cp -R /home/project1/Ecoli-Project1/Fastqs /your/destination/here/Fastqs
```

### 2. Executing the script
Once you have all of the needed files (and dependencies!) downloaded in your primary workspace (make sure everything is in the same directory), the script is easy to call.
```
python pipeline.py
```
The pipeline will create a new directory called "pipeline_files" and output everything it creates into that space.

### 3. Finding output files
All output files will be in a directory called "pipeline_files". Specific output files of interest can be found at the following locations:

TSV of all output metrics: /pipeline_files/full_output_data.tsv

Heatmap of comparisons: /pipeline_files/3_Metric_Heatmap.png

Bar chart of comparisons: /pipeline_files/barchart.png

MultiQC report (pre-trimming): /pipeline_files/MultiQC_Results/multiqc_report.html

MultiQC report (post-trimming): /pipeline_files/Trimmomatic_Results/MultiQC_Results/multiqc_report.html
