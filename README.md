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

## Tools
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

## Installing the Softwares
Make sure all of the tools are installed within the directory you will running everything. 

In case you have issues with MultiQC, MASH, and Trimmomatic, please refer to this for additional help - also located within their githubs tagged up above. 

Downloading MultiQC: 
```bash

git clone https://github.com/MultiQC/MultiQC.git

cd MultiQC

pip install .
```

Downloading MASH:
```bash

wget https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar

tar -xvf mash-Linux64-v2.3.tar

cp mash-Linux64-v2.3/mash .
```

Downloading Trimmomatic:
```bash
Micky add stuff here.
```
## What this Pipeline is Doing and the Output it Gives
This pipeline is comparing 12 _E. coli_ strains to their ancestor strain, from which diverged into 12 different populations. To compare these strains to one another, we used 3 different metrics, that have not been used to compare these specific strains in other studies using the same strains. These metrics consisted of ANI, MASH, and dDDH. 

Previously we downloaded all of the SRR numbers from NCBI, and then ran them through Prefetch and Fasterq-Dump. This resulted in FASTQ files that could be used for running this pipeline. 

### 1. Discovering Read Quality and Trimming
The pipeline first starts off with using FastQC to access the quality of each of the raw reads. MultiQC was then used to give a summarized singular file with all the reads output from FastQC. This makes it easier to read, rather than having a separate file for each sample. 

The reads were then trimmed using Trimmomatic. Then, these trimmed reads were used again in FastQC and MultiQC to see the differences in read quality compared to what it was before they were trimmed. 

### 2. Assembling Genomes to Use for Comparisons
Using the trimmed reads, SpAdes will create assembled genomes that can be used to compare strains to one another. This will output FASTA files that will be used in the following 3 metrics. 

### 3. 3 Comparison Metrics
After obtaining the FASTA files, ANI, MASH, and dDDH, will take each of the FASTA files and compare them to the ancestor strain. This will give a total of 12 comparisons (each most recent population strain compared to the ancestor strain). This will output a CSV file for each one of these metrics. 

### 4. Comparing Metrics
To compare the results from the 3 different metrics, a TSV file with all 3 metric information will output.

### 5. Visualization
Using this made TSV file, it will create a heatmap and a bar chart showing the different similarity calculations between each of the metrics.

## TO TEST RUN - Use this sample data:
This data involves shortened fastqs that were taken from the original data reads: [Sample Data]().

## How to Run the Script

### 1. 




