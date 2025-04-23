import csv
import subprocess
import os

main_dir = "/home/2025/Project1" #main directory - change here 
os.makedirs(main_dir, exist_ok=True) #make the directory if not there
os.chdir(main_dir) #move into the directory

csv_strains = os.path.join(main_dir, "sample_data.csv") #csv with strain info
#sra_prefetch_output = os.path.join(main_dir, "Prefetch_files") #where the sra files will go from prefetch output
#fastqs_faster_output = os.path.join(main_dir, "Fastqs") #where fasterq output is going to go (fastqs)

#os.makedirs(sra_prefetch_output, exist_ok=True) #make directory if not there
#os.makedirs(fastqs_faster_output, exist_ok=True) #same here - make directory if not there

'''
# Step 1: opening csv file with all the strain info and moving the SRR numbers to a list

#empty list to hold all the accession numbers
accessions = []
groups = {
    "Ara+1": [],
    "Ara+2": [],
    "Ara+3": [],
    "Ara+4": [],
    "Ara+5": [],
    "Ara+6": [],
    "Ara-1": [],
    "Ara-2": [],
    "Ara-3": [],
    "Ara-4": [],
    "Ara-5": [],
    "Ara-6": []
}
with open(csv_strains, mode = 'r') as f: #opening and indicating reading the csv file
    csv_dict = csv.DictReader(f) #creating a dictionary for each row where the headers are the keys
    for row in csv_dict: #for each line in the csv dictionary
        accession = row["Accession"].strip() #getting the accession numbers and stripping extra whitespace
        accessions.append(accession) #and adding it to assessions list
        group = row["\ufeffPopulation"].strip()
        if group in groups:
            groups[group].append(accession)
        #print(row)
    
#printing to see that the accession numbers were extracted correctly 
print(accessions)
print(groups)


# Step 2: get fastqs from NCBI + run Prefetch and Fasterq-dump
for accession in accessions: #for each accession in the accessions list
    #print(type(accession))
    #prefetch gets the SRA number from NCBI
    prefetch = ["prefetch", accession, "-O", sra_prefetch_output] #outputting into sra_prefetch_output folder
    #print(prefetch)
    subprocess.run(prefetch) #running the prefetch command
    #subprocess.run(prefetch, check = True) #running the prefetch command

    #path with sra file that was downloaded by prefetch and labeling with accession number
    sra_pathway = os.path.join(sra_prefetch_output, accession, f"{accession}.sra")
    #make directories for fasterq output - and labeling with accession number
    output_dir = os.path.join(fastqs_faster_output, accession)
    os.makedirs(output_dir, exist_ok=True) #making directory if not already made

    #running fasterq to make the fastqs for each sra number
    fasterq = ["fasterq-dump", sra_pathway, "-O", output_dir]
    #running the fasterq command
    subprocess.run(fasterq)
    #subprocess.run(fasterq, check = True)

print("All fastqs have been downloaded.")

'''


sample_reads = os.path.join(main_dir, "sample_reads")
accessions = os.listdir(sample_reads)

# 1. FastQC

#function to find fastq
def get_fastq(folder): #will have to put the path of the fastqs that you wanted to use to get
    #making an empty list to put fastq names
        fastqs = []
        for file in os.listdir(folder): #for each file within the folder
            if file.endswith(".fastq"): #if the file ends with .fastq
                fastqs.append(file) #add it to the fastqs list
        return fastqs

#function to run fastqc
def fastqc(accessions, fastqs_faster_output): 
    for name in accessions: #for each strain name in the accessions list 
        faster_output = os.path.join(fastqs_faster_output, name) #going through fasterq dump output to get the fastq name

        #going through all the files in the directory for the fastq outputs
        fastq_files = get_fastq(faster_output) #using paht of where the fastq names are at

        #creating the place where the fastqc output is going to go
        fastqc_output = os.path.join(faster_output, "FastQC_Results") #path will be to the place where fasterq dump output is, adding a folder for fastqc output
        os.makedirs(fastqc_output, exist_ok = True) #make the directory if not already there

        #now running each fastq on fastqc
        for fastq in fastq_files: #for each fastq in fastq_files list
            fastq_get = os.path.join(faster_output, fastq) #getting each of the fastqs by their name from the path variable faster_output
            #run fastqc 
            fastqc_command = ["fastqc", fastq_get, "-o", fastqc_output] #input the fastqs and output into the fastqc output place
            subprocess.run(fastqc_command, check = True) #running the fastqc command
    
    #so we know this step is done
    print("FastQC is done for ALL the files.")

#2. MultiQC
def multiqc(input_dir): #function to run multiqc
    
    #making the output directory for multiqc results
    multiqc_output = os.path.join(main_dir, "MultiQC_Results") #main directory pathway (Ecoli-Project1) and adding MultiQC_results section
    os.makedirs(multiqc_output, exist_ok=True) #making folder if not already there

    #multiqc command - will look through input that it gets and summarize all the fastqc results that it gets
    multiqc_command = ["multiqc", input_dir, "-o", multiqc_output]

    #running the multiqc command
    subprocess.run(multiqc_command, check = True)

    #making sure this part was done
    print(f"MultiQC is done :) .")

#Running fastqc and multiqc on the fastq files
fastqc(accessions, sample_reads) #calling fastqc to run with the sample reads output
multiqc(sample_reads) #calling multiqc to run wiht the fastqc output 

def run_trimmomatic(): # Define the directory containing input FASTQ files
    fastqs_dir = os.path.join(main_dir, "Fastqs")
    trim_output = os.path.join(main_dir, "Trimmomatic_Results")
    os.makedirs(trim_output, exist_ok=True) # Define the output directory where trimmed FASTQ files will be saved

  # Define paths using main_dir so it's relative to the project directory
    trimmomatic_dir = os.path.join(main_dir, "..", "miniconda3", "envs", "trimmomatic_env", "share", "trimmomatic-0.39-2")
    trimmomatic_jar = os.path.join(trimmomatic_dir, "trimmomatic.jar")
    adapter_pe = os.path.join(trimmomatic_dir, "adapters", "TruSeq3-PE.fa")
    adapter_se = os.path.join(trimmomatic_dir, "adapters", "TruSeq3-SE.fa")
    
    # Loop through each accession directory inside the sample_reads folder
    for accession in os.listdir(fastqs_dir):
        accession_dir = os.path.join(fastqs_dir, accession)
        if not os.path.isdir(accession_dir):
            continue # Skip if it's not a directory

        files = os.listdir(accession_dir)  # Try to identify paired-end FASTQ files (R1 and R2)
        r1 = next((f for f in files if "_R1" in f and f.endswith(".fastq")), None)
        r2 = next((f for f in files if "_R2" in f and f.endswith(".fastq")), None)
    # If both R1 and R2 are found, run Trimmomatic in paired-end (PE) mode
        if r1 and r2:
            input_r1 = os.path.join(accession_dir, r1)
            input_r2 = os.path.join(accession_dir, r2)
            cmd = [
                "java", "-jar", trimmomatic_jar, "PE", "-phred33",
                input_r1, input_r2,
                os.path.join(trim_output, f"{accession}_R1_paired.fastq"),
                os.path.join(trim_output, f"{accession}_R1_unpaired.fastq"),
                os.path.join(trim_output, f"{accession}_R2_paired.fastq"),
                os.path.join(trim_output, f"{accession}_R2_unpaired.fastq"),
                f"ILLUMINACLIP:{adapter_pe}:2:30:10", "LEADING:3",
                "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"]
             # Adapter trimming settings, Trim leading low-quality bases (below quality 3), Trim trailing low-quality bases (below quality 3)
             # Sliding window trimming with window size 4, quality threshold 15
             # Drop reads shorter than 36 bases
            print(f" Running PE Trimmomatic on {accession}")
        else:
             # If paired-end not detected, default to single-end (SE) mode
            se_file = next((f for f in files if f.endswith(".fastq")), None)
            if not se_file:
                print(f"No FASTQ found for {accession}")
                continue
            input_se = os.path.join(accession_dir, se_file)
            cmd = [
                "java", "-jar", trimmomatic_jar, "SE", "-phred33",
                input_se, os.path.join(trim_output, f"{accession}_trimmed.fastq"),
                f"ILLUMINACLIP:{adapter_se}:2:30:10", "LEADING:3",
                "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"
            ]
            print(f" Running SE Trimmomatic on {accession}")
 # run the constructed Trimmomatic command
        subprocess.run(cmd, check=True)
        print(f" Finished Trimmomatic for {accession}")

    print("Trimmomatic is done.")
    return trim_output

trim_output = run_trimmomatic()