# Running this Pipeline

## Dependencies
Python modules:
* csv
* subprocess
* os
* more TBD

Tools:
* [FastQC](https://github.com/s-andrews/FastQC)
* [Trimmomatic](https://github.com/timflutre/trimmomatic)
* [SPAdes](https://github.com/ablab/spades)
* [FastANI](https://github.com/ParBLiSS/FastANI)
* [MASH](https://github.com/marbl/mash)
## Data Retrieval 
* Download raw sequencing reads of 264 LTEE genomes from NCBI SRA
* Retrieve ancestral strain genome [REL606](https://www.ncbi.nlm.nih.gov/nuccore/NC_012967.1)

## Input Data:
* 256 _E. coli_ genomes from the Long Term Evolution Experiment (LTEE)
* One ancestor (REL606) genome from the LTEE
## References 
* [Wiki Link](https://github.com/lexrex333/Ecoli-Project1/wiki) for project 
* [Introduction to the Long-Term Evolution Experiment (LTEE)](https://the-ltee.org/about/) 
* [Tempo and mode of genome evolution in a 50,000-generation experiment](https://www.nature.com/articles/nature18959)

## Installing the Softwares
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
