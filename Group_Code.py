import csv
import subprocess
import os
import pandas as pd
from itertools import combinations


main_dir = "/home/project1/Ecoli-Project1" #main directory - change here 
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

#3. Trimmomatic

def run_trimmomatic():
    fastqs_dir= os.path.join(main_dir, "sample_reads") #set directory containing raw fastq files  
    trim_output= os.path.join(main_dir, "Trimmomatic_Results") #set directory to store the trimmed output files
   
    os.makedirs(trim_output, exist_ok=True) #make sure the output directory exsists                                    
    # Define the paths to the trimmomatic jar and adapter file 
    trimmomatic_jar = os.path.join(main_dir, "Trimmomatic-0.39", "trimmomatic-0.39.jar")
    adapter_file = os.path.join(main_dir, "Trimmomatic-0.39", "adapters", "TruSeq3-SE.fa")
    
    for root, dirs, files in os.walk(fastqs_dir): #go through the fastq directory and subdirectories to find file 
     for file in files:  #loop through all files ends with .fastq
        if file.endswith(".fastq"):
            fastq_path = os.path.join(root, file)  #create full path to the file
            #construct the path for the trimmed output file
            trimmed_fastq = os.path.join(trim_output, f"{file}")

            #construct the Trimmomatic command for trimming the fastq files 
            #run trimmomatic, single end reads, specify phred quality score, input Fastq file path and output trimmed file path
            #remove leading and trailing low quality bases, trimming if average quality <15, only retain reads that are 36 bases or longer 
            trimmomatic_command = [
                "/home/project1/miniconda3/envs/trimmomatic_env/bin/trimmomatic" , "SE", "-phred33", fastq_path, trimmed_fastq,
                 f"ILLUMINACLIP:{adapter_file}:2:30:10", "LEADING:3", "TRAILING:3",
                "SLIDINGWINDOW:4:15", "MINLEN:36"
            ] 
            print("Running:", " ".join(trimmomatic_command))
 #run the Trimmomatic command
            subprocess.run(trimmomatic_command, check=True)

     #print a message indicating that Trimmomatic processing is complete for all files
    print("Trimmomatic is done.")
    return trim_output

trim_output = run_trimmomatic() #running trimmomatic

#another function to run fastqc because it didn't want to work for trimmed files :(
def fastqc_trimmed(folder): 
    fastq_files = get_fastq(folder) #calling function to make a fastqs list

    #creating the output folder
    fastqc_output = os.path.join(folder, "FastQC_Results") #where the results will go -- within the trimmomatic results
    os.makedirs(fastqc_output, exist_ok=True) #making directory if not already there

    for fastq in fastq_files: #for each of the fastqs within the fastq files list
        fastq_path = os.path.join(folder, fastq) #path of where the fastqs are at 
        fastqc_command = ["fastqc", fastq_path, "-o", fastqc_output] #fastqc command
        subprocess.run(fastqc_command, check = True) #run fastqc command 

    
#4. FastQC again
#have to run FastQC and MultiQC on the trimmed files 
fastqc_trimmed(trim_output) #calling fastqc to run with the trimmed output
multiqc(trim_output) #calling multiqc to run wiht the fastqc output just made from the trimmed fastqs



#5. SPAdes

read_types = {} #create dictionary to hold read types (single/paired) for each genome

#os.chdir("Ecoli-Project1") #navigate to folder containing csv file

with open("sample_data.csv") as file: #access data from csv file
    rd = csv.reader(file, delimiter=",")
    next(file) #skip first line
    for row in rd:
        read_types[f"{row[6]}"] = f"{row[7]}" #add SRR numbers as keys and "single"/"paired" as values in dictionary

#os.mkdir("sample_assemblies") #uncomment this out for the official one

#function to run spades given a dictionary of read numbers and single/paired read statuses
def run_spades(read_dict):
    for key, value in read_dict.items():
        # os.chdir(f"{key}") #navigate to folder containing the specific genome's read file(s)
        if value == "single":
            filename = f"Trimmomatic_Results/{key}"+".fastq" #name of file containing reads for this genome
            spades_command = f"spades.py -s {filename} -o sample_assemblies/{key}_assembly" #building command to run spades and output in the same folder where the reads are contained
        elif value == "paired":
            file1 = f"Trimmomatic_Results/{key}"+"_1.fastq" #names of files containing reads for this genome
            file2 = f"Trimmomatic_Results/{key}"+"_2.fastq"
            spades_command = f"spades.py -1 {file1} -2 {file2} -o sample_assemblies/{key}_assembly" #building command to run spades
        os.system(spades_command) #execute command
        
print(run_spades(read_types))
print("SPAdes complete")
'''

#6. FastANI
'''
#function to keep running if fastas are too short
def keep_running(fasta_path, min_len = 50):
    with open(fasta_path, 'r') as f: #opening in read mode
        seq = '' #holding sequence
        for line in f: #for each line in the file
            line = line.strip() #stripping whitespace
            if line.startswith('>'): #if line starts with >
                if len(seq) >= min_len: #if the length is greater than or equal to min length
                    return True #return as true
                seq = ''  #resetting for the next contig
            else: #if its a sequence line,  
                seq += line #add it to the current sequence
        return len(seq) >= min_len #checking the last contig
    '''
#this function is going to run fastani within each pop group and compares each sample to the ancestor strain
def fastani_run(main_dir, sample_csv, fastani_output="FastANI_Results", ancestor_name = "REL606A.fasta"):

    #place to save fastani output
    fastani_output = os.path.join(main_dir, 'FastANI_Results')
    #make folder if not already made
    os.makedirs(fastani_output, exist_ok=True)

    #loading sample data to look off of
    df = pd.read_csv(sample_csv)
    #debugging based on how it reads it
    df.columns = df.columns.str.strip()
    #debug
    df["Accession"] = df["Accession"].str.strip()
    df["Population"] = df["Population"].str.strip()

    sample_fastas = {} #place to put paths of each spades fasta file

    #dictionary to hold the population groups
    pop_groups = {}

    #place to put all FastANI results before putting it to csv file
    results = [] #empty list

    #going through each row of the csv to group the fastas based on their groups
    for _, row in df.iterrows(): #_ ignores the index and just gets the rows
        sample = row["Accession"].strip() #getting the sample name
        population = row["Population"].strip() #getting population row
        
        #path to spades fastas
        fasta_path = os.path.join(main_dir,"shitty_sample_assemblies", f"{sample}_assembly", "contigs.fasta") #fasta path to spades output file
        
        if os.path.isfile(fasta_path): #if there actually is a file in that path
            sample_fastas[sample] = fasta_path #add the folder name and path to the fasta file to the fastas dictionary 

            #adding the sample to the pop group that it belongs to
            if population not in pop_groups: #if the population is already not added to the groups dictionary
                pop_groups[population] = [] #add the population as a key
            pop_groups[population].append(sample) #add the sample to the population 
        
        #debugging 
        else:
            print(f"No fasta for {sample} at {fasta_path}")

    #debugging
    print(list(sample_fastas.keys()))
    print(pop_groups.keys())

    #Step ONE: comparing within each population
    for population, samples in pop_groups.items():

        #getting all the comparisons bw 2 samples between this individual group
        #combinations(samples,2) is going to do all unique pairs without repeating them or comparing them to themselves
        for sample1, sample2 in combinations(samples, 2):
            fasta1 = sample_fastas[sample1]
            fasta2 = sample_fastas[sample2]

            #making sure it still runs with the fastas we have
           # if not keep_running(fasta1) or not keep_running(fasta2):
               # print(f"Skipping {sample1} vs {sample2}- not enough contigs - too short.")
                #continue

            #making a temp file for fastani output
            fastani_temp = os.path.join(fastani_output, f"{sample1}_vs_{sample2}.txt")

            #fastani command -- will need to change the pathway by manually getting fastANI
            fastani = [os.path.join(main_dir, "fastANI"), "-q", fasta1, "-r", fasta2, "--fragLen", "500", "-o", fastani_temp] #storing output to terminal instead to write to csv
            #running the command
            subprocess.run(fastani, check=True)
            
            #debugging - to get proof that there isn't that much of a match
            if os.path.getsize(fastani_temp) == 0:
                print(f"No ANI result written for {sample1} vs {sample2}.")
            else: 
            
                #going through the output and adding to the results list
                with open(fastani_temp, 'r') as f:
                    for line in f:
                        _, _, ani, fragments, total = line.strip().split("\t") #strip whitespace and split by tab
                        results.append({
                            "Query": sample1,
                            "Reference": sample2,
                            "ANI": ani,
                            "Matching_Fragments": fragments,
                            "Total_Fragments": total,
                            "Group": population,
                            "Comparison_Method": "Within Population"
                        })

    #STEP TWO: Comparing ALL the samples to the ancestor
    ancestor = sample_fastas.get(ancestor_name) #might have to change this later when added to pipeline

    if not ancestor: #if ancestor name not found in the sample list
        #debugging
        print(f"Ancestor {ancestor_name} not found.")
        return
    
    #comparing each fasta to the ancestor except the ancestor fasta
    for sample, fasta in sample_fastas.items(): #going through sample_fastas dictionary
        if sample == ancestor_name: #if the sample equals to the ancestor
            continue #ignore and keep going

        #to make sure it keeps running even if some are empty
       # if not keep_running(fasta) or not keep_running(ancestor):
          #  print(f"Skipping {sample} vs Ancestor - fasta too short.")
          #  continue

        #place to put fastani output temporarily
        fastani_temp2 = os.path.join(fastani_output, f"{sample}_vs_{ancestor_name}.txt")
        #the fastani command
        fastani_command = [os.path.join(main_dir, "fastANI"), "-q", fasta, "-r", ancestor, "--fragLen", "500", "-o", fastani_temp2]
        #running the fastani command
        subprocess.run(fastani_command, check = True)

        #this looks through the csv and gets the row where the accession is equal to the sample
        #and gets the value that is in the population column that matches that 
        #in short talk, it is just getting the pop group that the sample belongs to
        group = df[df["Accession"] == sample]["Population"].values[0]

        #debugging - to get proof that there isn't that much of a match
        if os.path.getsize(fastani_temp) == 0:
            print(f"No ANI result written for {sample1} vs {sample2}.")
        else: 

            #going through the output and adding to the results list
            with open(fastani_temp2, 'r') as f:
                for line in f:
                    _, _, ani, fragments, total = line.strip().split("\t") #strip whitespace and split by tab
                    results.append({"Query": sample,
                                    "Reference": ancestor_name,
                                    "ANI": ani,
                                    "Matching_Fragments": fragments,
                                    "Total_Fragments": total,
                                    "Group": group,
                                    "Comparison_Method": "To Ancestor"
                                    })
    
    results_to_df = pd.DataFrame(results) #changing the results to a dataframe
    results_to_df.to_csv(os.path.join(fastani_output, "FastANI_Results_ALL.csv"), index = False)
    
    print("FastANI has been completed!!")

fastani_run(main_dir, sample_csv = os.path.join(main_dir, "sample_data.csv"), fastani_output = "FastANI_Results", 
            ancestor_name = "SRR22764941")

