import csv
import subprocess
import os
import pandas as pd
from itertools import combinations
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns

#create directory to hold ALL FILES that this pipeline will create and work with
os.mkdir("pipeline_files")
os.system("cp recent_strains.csv pipeline_files/recent_strains.csv")
os.system("cp -R Fastqs pipeline_files/Fastqs") #optional: only if using sample data
os.chdir("pipeline_files")


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
with open("recent_strains.csv", mode = 'r') as f: #opening and indicating reading the csv file
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

#set main directory as current working directory
main_dir = os.getcwd()
fastqs = os.path.join(main_dir, "Fastqs")
print(fastqs)
accessions = os.listdir(fastqs)
print(accessions)


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

        #now running each fastq on fastqc
        for fastq in fastq_files: #for each fastq in fastq_files list
            fastq_get = os.path.join(faster_output, fastq) #getting each of the fastqs by their name from the path variable faster_output
            #run fastqc 
            fastqc_command = ["fastqc", fastq_get, "-o", "FastQC_Results"] #input the fastqs and output into the fastqc output place
            subprocess.run(fastqc_command, check = True) #running the fastqc command
    
    #so we know this step is done
    print("FastQC is done for ALL the files.")

#2. MultiQC

QCresults = os.path.join(main_dir, "FastQC_Results") #read results from FastQC folder

def multiqc(input_dir): #function to run multiqc
    
    #making the output directory for multiqc results
    os.mkdir("MultiQC_Results") #main directory pathway and adding MultiQC_results section

    #multiqc command - will look through input that it gets and summarize all the fastqc results that it gets
    multiqc_command = ["multiqc", input_dir, "-o", "MultiQC_Results"]

    #running the multiqc command
    subprocess.run(multiqc_command, check = True)

    #making sure this part was done
    print(f"MultiQC is done :) .")

# output folder for results
os.mkdir("FastQC_Results")
# Running fastqc and multiqc on the fastq files
fastqc(accessions, fastqs) #calling fastqc to run with the sample reads output
multiqc(QCresults) #calling multiqc to run with the FastQC output 



# if using sample reads, skip over trimming and second round of QC
# Trimmomatic would trim short reads beyond usability and SPAdes would not be able to assemble

with open("Fastqs/SRR2584406/SRR2584406_1.fastq") as f:
    x = len(f.readlines()) #save length of readfile as x

#proceed through this section only if reads are longer than sample length
if x > 120010:
    #3. Trimmomatic

    def run_trimmomatic():
        # Define the directory containing input FASTQ files
        fastqs_dir = os.path.join(main_dir, "Fastqs")
        os.mkdir("Trimmomatic_Results")
        trim_output = os.path.join(main_dir, "Trimmomatic_Results")
        os.makedirs(trim_output, exist_ok=True)  # Define the output directory where trimmed FASTQ files will be saved

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
            r1 = next((f for f in files if "_1" in f and f.endswith(".fastq")), None)
            r2 = next((f for f in files if "_2" in f and f.endswith(".fastq")), None)
        # If both R1 and R2 are found, run Trimmomatic in paired-end (PE) mode
            if r1 and r2:
                input_r1 = os.path.join(accession_dir, r1)
                input_r2 = os.path.join(accession_dir, r2)
                cmd = [
                    "java", "-jar", trimmomatic_jar, "PE", "-phred33",
                    input_r1, input_r2,
                    os.path.join(trim_output, f"{accession}_1_paired.fastq"),
                    os.path.join(trim_output, f"{accession}_1_unpaired.fastq"),
                    os.path.join(trim_output, f"{accession}_2_paired.fastq"),
                    os.path.join(trim_output, f"{accession}_2_unpaired.fastq"),
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
                    input_se, os.path.join(trim_output, f"{accession}.fastq"),
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



    #another function to run fastqc because it didn't want to work for trimmed files :(
    def fastqc_trimmed(folder): 
        fastq_files = get_fastq(folder) #calling function to make a fastqs list

        #creating the output folder
        os.chdir("Trimmomatic_Results")
        os.mkdir("FastQC_Results")#where the results will go -- within the trimmomatic results

        for fastq in fastq_files: #for each of the fastqs within the fastq files list
            fastq_path = os.path.join(folder, fastq) #path of where the fastqs are at 
            fastqc_command = ["fastqc", fastq_path, "-o", "FastQC_Results"] #fastqc command
            subprocess.run(fastqc_command, check = True) #run fastqc command 

        
    #4. FastQC again
    #have to run FastQC and MultiQC on the trimmed files 
    fastqc_trimmed(trim_output) #calling fastqc to run with the trimmed output
    multiqc(trim_output) #calling multiqc to run wiht the fastqc output just made from the trimmed fastqs

    #move back to main directory
    os.chdir("..")

#5. SPAdes

read_types = {} #create dictionary to hold read types (single/paired) for each genome

with open("recent_strains.csv") as file: #access data from csv file
    rd = csv.reader(file, delimiter=",")
    next(file) #skip first line
    for row in rd:
        read_types[f"{row[6]}"] = f"{row[7]}" #add SRR numbers as keys and "single"/"paired" as values in dictionary

os.mkdir("assemblies")

#function to run spades given a dictionary of read numbers and single/paired read statuses
def run_spades(read_dict):
    for key, value in read_dict.items():
        if x > 120010: #if reads are full-length, pull from trimmomatic results
            if value == "single":
                filename = f"Trimmomatic_Results/{key}"+".fastq" #name of file containing reads for this genome
                spades_command = f"spades.py -s {filename} -o assemblies/{key}_assembly" #building command to run spades and output in the same folder where the reads are contained
            elif value == "paired":
                file1 = f"Trimmomatic_Results/{key}"+"_1_paired.fastq" #names of files containing reads for this genome
                file2 = f"Trimmomatic_Results/{key}"+"_2_paired.fastq"
                spades_command = f"spades.py -1 {file1} -2 {file2} -o assemblies/{key}_assembly" #building command to run spades
        elif x < 120010: #if reads are shortened, pull from fastqs
            if value == "single":
                filename = f"Fastqs/{key}/{key}"+".fastq" #name of file containing reads for this genome
                spades_command = f"spades.py -s {filename} -o assemblies/{key}_assembly" #building command to run spades and output in the same folder where the reads are contained
            elif value == "paired":
                file1 = f"Fastqs/{key}/{key}"+"_1.fastq" #names of files containing reads for this genome
                file2 = f"Fastqs/{key}/{key}"+"_2.fastq"
                spades_command = f"spades.py -1 {file1} -2 {file2} -o assemblies/{key}_assembly" #building command to run spades            
        os.system(spades_command) #execute command
        
print(run_spades(read_types))
print("SPAdes complete")


#6. FastANI
#this function is going to run fastani within each pop group (even tho we didnt do this
# but if wanted to use more than one sample in each population can do so)
# and compares each sample to the ancestor strain
def fastani_run(main_dir, sample_csv, fastani_output="FastANI_Results", ancestor_name = "REL606A.fasta"):

    #place to save fastani output
    fastani_output = os.path.join(main_dir, 'FastANI_Results')
    #make folder if not already made
    os.makedirs(fastani_output, exist_ok=True)

    #loading sample data to look off of
    df = pd.read_csv(sample_csv)
    #debugging based on how it reads it
    df.columns = df.columns.str.strip()
    #debug - taking out whitespace and putting string
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
        fasta_path = os.path.join(main_dir,"assemblies", f"{sample}_assembly", "contigs.fasta") #fasta path to spades output file
        
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

            #making a temp file for fastani output
            fastani_temp = os.path.join(fastani_output, f"{sample1}_vs_{sample2}.txt")

            #fastani command -- will need to change the pathway by manually getting fastANI
            fastani = [os.path.join(main_dir, "..", "fastANI"), "-q", fasta1, "-r", fasta2, "--fragLen", "1000", "-o", fastani_temp] #storing output to terminal instead to write to csv
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

        #place to put fastani output temporarily
        fastani_temp2 = os.path.join(fastani_output, f"{sample}_vs_{ancestor_name}.txt")
        #the fastani command
        fastani_command = [os.path.join(main_dir, "..", "fastANI"), "-q", fasta, "-r", ancestor, "--fragLen", "1000", "-o", fastani_temp2]
        #running the fastani command
        subprocess.run(fastani_command, check = True)

        #this looks through the csv and gets the row where the accession is equal to the sample
        #and gets the value that is in the population column that matches that 
        #in short talk, it is just getting the pop group that the sample belongs to
        group = df[df["Accession"] == sample]["Population"].values[0]

        #debugging - to get proof that there isn't that much of a match
        if os.path.getsize(fastani_temp2) == 0:
            print(f"No ANI result written.")
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

#calling the function to run
fastani_run(main_dir, sample_csv = os.path.join(main_dir, "recent_strains.csv"), fastani_output = "FastANI_Results", ancestor_name = "SRR22764941")


#7. MASH
Sketch_Output = "Sketches_for_MASH"
if os.path.exists(Sketch_Output):  # Force fresh sketch with updated paths
        os.rmdir(Sketch_Output)
os.mkdir(Sketch_Output)
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
    "Ara-6": [],
}

homedir = os.getcwd()
base_path = os.path.join(homedir, "assemblies")
results_file = "mash_results.tsv"
mash_path = os.path.join(homedir, "..", "mash")

# Clearing output file so that each run produces a new file output
with open(results_file, "w") as f:
    f.write("Sample_1\tSample_2\tMASH_Distance\tSimilarity\tpvalue\tShared_Hashes\n")

# Reading in Data
with open("recent_strains.csv", mode='r') as f:
    csv_dict = csv.DictReader(f)
    for row in csv_dict:
        accession = row["Accession"].strip()
        accessions.append(accession)
        group = row["\ufeffPopulation"].strip()
        if group in groups:
            groups[group].append(accession)

# Add ancestor accession to each group for comparison
ancestor = "SRR22764941"
for group in groups:
    groups[group].append(ancestor)

# Get FASTQ file path for analysis
def get_contigs_path(accession):
    return os.path.join(base_path, f"{accession}_assembly", "contigs.fasta"), accession

# Create a sketch for an accession, return the .msh file path, .msh will be stored with SPADES output
def create_sketch(contigs_path,SRR):
    sketch_path = Sketch_Output + "/" + SRR + ".msh"
    if os.path.exists(sketch_path):  # Force fresh sketch with updated paths
        os.remove(sketch_path)
    #if not os.path.exists(sketch_path):  # Only sketch if it doesn't exist
    sketch_command = [mash_path, "sketch", "-o", sketch_path, contigs_path]
    subprocess.run(sketch_command, check=True)
    return sketch_path

# Run MASH and write results to file
def MASH(one_accession, two_accession):
    one_contigs, one_srr = get_contigs_path(one_accession)
    two_contigs, two_srr = get_contigs_path(two_accession)

    sketch_one = create_sketch(one_contigs, one_srr)
    sketch_two = create_sketch(two_contigs, two_srr)

    mash_command = [mash_path, "dist", sketch_one, sketch_two]
    result = subprocess.run(mash_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

  # Parse the MASH output
    if result.stdout:
        parts = result.stdout.strip().split("\t")
        if len(parts) == 5:
            distance = parts[2]
            similarity = 1 - float(distance)
            pvalue = parts[3]
            shared_hashes = parts[4]

            with open(results_file, "a") as f:
                f.write(f"{one_accession}\t{two_accession}\t{distance}\t{similarity}\t{pvalue}\t{shared_hashes}\n")
        else:
            print(f"Unexpected MASH output format for {one_accession} vs {two_accession}: {result.stdout}")
    else:
        print(f"Error: No output from MASH for {one_accession} vs {two_accession}")
        
# Loop through all pairwise comparisons
for group_name, accession_list in groups.items():
    for one, two in combinations(accession_list, 2):
        print(one,two)
        MASH(one, two)


#8. dDDH

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

with open("recent_strains.csv", mode = 'r') as f: #opening and reading the csv file
    csv_dict = csv.DictReader(f) #creating a dictionary for each row where the headers are the keys
    for row in csv_dict: #for each line in the csv dictionary
        accession = row["Accession"].strip() #getting the accession numbers and stripping extra whitespace
        group = row["\ufeffPopulation"].strip()
        if group in groups:
            groups[group].append(accession)
    for group in groups:
        groups[group].append("SRR22764941") #add ancestor accession to the end of each population group
print(groups)


def run_blast(comparedict):
    os.mkdir("blast_output")
    #create a tsv file to hold output data that we can work with
    f = open("blast_output_data.tsv", "w")
    f.write("sample 1\tsample 2\tdDDH\tpopulation\n")
    f.close()
    for key, value in comparedict.items():
        for i in range(len(value)):
            if i+1 == len(value):
                continue
            for j in range(1, len(value)):
                if i+j >= len(value):
                    continue
                print(value[i], value[i+j])
                blast_command = f"blastn -query assemblies/{value[i]}_assembly/contigs.fasta -subject assemblies/{value[i+j]}_assembly/contigs.fasta -max_hsps 1 -out blast_output/{value[i]}_{value[i+j]}_blast.tsv -outfmt '6 nident length'"
                os.system(blast_command) #execute blast command
                with open(f"blast_output/{value[i]}_{value[i+j]}_blast.tsv") as file: #access blast results
                    rd = csv.reader(file, delimiter="\t")
                    lengthlist = []
                    nidentlist = []
                    for row in rd: #pull nident & length values for each contig alignment from blast results
                        nident = int(row[0])
                        nidentlist.append(nident)
                        length = int(row[1])
                        lengthlist.append(length)
                total_nident = 0
                total_length = 0
                for x in range(len(nidentlist)): #add up all nident and length values
                    total_nident += nidentlist[x]
                    total_length += lengthlist[x]
                dDDH = (total_nident / total_length) #divide total nident / total length to get overall dDDH
                print(dDDH)
                f = open("blast_output_data.tsv", "a")
                f.write(f"{value[i]}\t{value[i+j]}\t{dDDH}\t{key}\n")
                f.close()
    return("dDDH complete")

print(run_blast(groups))



#9. Putting the CSVs into one TSV

print("Creating results file...")

f = open("full_output_data.tsv", "w") #create an output file
f.close()

#make lists to hold all the important imformation from these files
sample1list = []
sample2list = []
ANIlist = []
MASHlist = []
MASHsimlist = []
dDDHlist = []
poplist = []

#read dDDH file and pull information
with open("blast_output_data.tsv") as file: #access data from tsv file
    rd = csv.reader(file, delimiter="\t")
    next(file) #skip first line
    for row in rd:
        sample1list.append(f"{row[0]}")
        sample2list.append(f"{row[1]}")
        dDDHlist.append(f"{row[2]}")
        poplist.append(f"{row[3]}")

#read MASH file and pull information
with open("mash_results.tsv") as file:
    rd = csv.reader(file, delimiter="\t")
    next(file) #skip first line
    for row in rd:
        MASHlist.append(f"{row[2]}")
        MASHsimlist.append(f"{row[3]}")

#read ANI file and pull information
with open("FastANI_Results/FastANI_Results_ALL.csv") as file:
    rd = csv.reader(file, delimiter=",")
    next(file) #skip first line
    for row in rd:
        ANIlist.append(f"{row[2]}")

#write info in new output file
f = open("full_output_data.tsv", "a")
f.write("sample_1\tsample_2\tANI\tMASH_dist\tMASH_sim\tdDDH\tpopulation\n") #create header row
for i in range(len(ANIlist)): #put info for each pair in each row
    f.write(f"{sample1list[i]}\t{sample2list[i]}\t{ANIlist[i]}\t{MASHlist[i]}\t{MASHsimlist[i]}\t{dDDHlist[i]}\t{poplist[i]}\n")
f.close()

print("Results file created")


# 10. Now Making Visuals
# 10a. Heatmap
print("Creating heatmap...")
#Step 1. Read the tsv file of the 3 metrics 
ani = pd.read_csv(os.path.join(main_dir, "full_output_data.tsv"), sep="\t") #tab separated

#renaming the columns to match what I got going on
ani = ani.rename(columns={
    "sample_1": "Query", #the first strain
    "sample_2": "Reference", #the strain it is being compared to - ancestor
    "population": "Group" #the population it is comign from 
})

#Step 2. CLeaning everything
#Keeping only the comments I need and want in the heatmap
summary = ani[["Group", "ANI", "MASH_sim", "dDDH"]].copy()

#grouping by population 
pop_summary = summary.groupby("Group").mean()

#multiply MASH and dDDH to match the ani numbers 
pop_summary["MASH_sim"] = pop_summary["MASH_sim"] * 100 
pop_summary["dDDH"] = pop_summary["dDDH"] * 100

#Putting the populations in order - alphabetically
pop_summary = pop_summary.sort_index()

#Step 3. PLotting the heatmap
plt.figure(figsize=(8, 10)) #figure - size  

#plotting the heatmap 
sns.heatmap( #using seaborn 
    pop_summary,
    cmap="viridis", #the colors for the heatmap
    annot=True, #showing the number inside each square 
    fmt=".3f", #3 decimal point for the numbers 
    cbar_kws={"label": "Similarity Percentage"} #label for colorbar
)

#titles and axis labels 
plt.title("ANI, MASH Similarity, and dDDH Compared to Ancestor", fontsize=16) 
plt.ylabel("Population", fontsize=14) #y-axis label
plt.xlabel("Metric", fontsize=14) #x-axis label

#making the layout clean
plt.tight_layout()

#saving the figure 
out_path = os.path.join(main_dir, "3_Metric_Heatmap.png")
plt.savefig(out_path, dpi=300) #saving the figure 
plt.close() #closing

print("Heatmap created")


#10b. Bar plot

print("Creating bar chart...")

#set size of plot
barWidth = 0.25
fig = plt.subplots(figsize =(12, 8)) 

#lists to hold data
ANI = []
MASH_1 = []
dDDH_1 = []

#pull data from tsv
with open("full_output_data.tsv") as file:
    rd = csv.reader(file, delimiter="\t")
    next(file) #skip first line
    for row in rd:
        ANI.append(float(f"{row[2]}")) #save data in lists
        MASH_1.append(float(f"{row[4]}"))
        dDDH_1.append(float(f"{row[5]}"))
MASHlist = []
dDDH = []
#convert MASH and dDDH decimals to percentages
for item in MASH_1:
    percent = float(item) * 100
    MASHlist.append(float(percent))
for item in dDDH_1:
    percent = float(item) * 100
    dDDH.append(float(percent))

#set bar sizes
br1 = np.arange(len(ANI)) 
br2 = [x + barWidth for x in br1] 
br3 = [x + barWidth for x in br2] 

#set colors and datasets for each bar
plt.bar(br1, ANI, color ='lightcoral', width = barWidth, 
        edgecolor ='grey', label ='ANI') 
plt.bar(br2, MASHlist, color ='khaki', width = barWidth, 
        edgecolor ='grey', label ='MASH') 
plt.bar(br3, dDDH, color ='skyblue', width = barWidth, 
        edgecolor ='grey', label ='dDDH') 

#set limits of y axis to zoom in on higher percentages
a, b = 75, 101
plt.ylim((a, b))

#set labels for titles, axes, and ticks
plt.title("Recent Population Members' Similarities to Ancestor REL606", fontweight ='bold', fontsize = 16)
plt.xlabel('Population', fontweight ='bold', fontsize = 12) 
plt.ylabel('Score (%)', fontweight ='bold', fontsize = 12) 
plt.xticks([r + barWidth for r in range(len(ANI))], 
        ['Ara+1', 'Ara+2', 'Ara+3', 'Ara+4', 'Ara+5', 'Ara+6', 'Ara-1', 'Ara-2', 'Ara-3', 'Ara-4', 'Ara-5', 'Ara-6'])
plt.legend(loc="lower right") #put legend box in bottom right
#save plot as png
plt.savefig('barchart.png')

print("Bar chart completed")

print("All output files saved")

print("Pipeline is complete!!!")