import csv
import subprocess
import os



main_dir = "/home/project1/Ecoli-Project1" #main directory - change here 
os.makedirs(main_dir, exist_ok=True) #make the directory if not there
os.chdir(main_dir) #move into the directory

csv_strains = os.path.join(main_dir, "LTEE_strains.csv") #csv with strain info
sra_prefetch_output = os.path.join(main_dir, "Prefetch_files") #where the sra files will go from prefetch output
fastqs_faster_output = os.path.join(main_dir, "Fastqs") #where fasterq output is going to go (fastqs)

os.makedirs(sra_prefetch_output, exist_ok=True) #make directory if not there
os.makedirs(fastqs_faster_output, exist_ok=True) #same here - make directory if not there


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


#3. Trimmomatic

def run_trimmomatic():
    fastqs_dir= os.path.join(main_dir, "Fastqs")
    trim_output= os.path.join(main_dir, "Trimmomatic_Results")
   
    os.makedirs(trim_output, exist_ok=True) #make sure the output directory exsists                                    
    for file in os.listdir(fastqs_dir): #loop through each file in the fastqs directory
        if file.endswith(".fastq"):
            #construct the full path to the fastq file
            fastq_path = os.path.join(fastqs_dir, file)
            #construct the path for the trimmed output file
            trimmed_fastq = os.path.join(trim_output, f"trimmed_{file}")

            #construct the Trimmomatic command
            trimmomatic_command = [
                "trimmomatic", "SE", "-phred33", fastq_path, trimmed_fastq,
                "ILLUMINACLIP:TruSeq3-SE.fa:2:30:10", "LEADING:3", "TRAILING:3",
                "SLIDINGWINDOW:4:15", "MINLEN:36"
            ]
 #run the Trimmomatic command
            subprocess.run(trimmomatic_command, check=True)

     #print a message indicating that Trimmomatic processing is complete for all files
    print("Trimmomatic is done.")
    return trim_output

trim_output = run_trimmomatic() #running trimmomatic

#4. FastQC again
#have to run FastQC and MultiQC on the trimmed files 
fastqc(accessions, trim_output) #calling fastqc to run with the trimmed output
multiqc(trim_output) #calling multiqc to run wiht the fastqc output just made from the trimmed fastqs



#5. SPAdes

read_types = {} #create dictionary to hold read types (single/paired) for each genome

os.chdir("Ecoli-Project1") #navigate to folder containing csv file

with open("LTEE_strains.csv") as file: #access data from csv file
    rd = csv.reader(file, delimiter=",")
    next(file) #skip first line
    for row in rd:
        read_types[f"{row[6]}"] = f"{row[7]}" #add SRR numbers as keys and "single"/"paired" as values in dictionary

os.chdir("Fastqs") #navigate to folder containing genomes

#function to run spades given a dictionary of read numbers and single/paired read statuses
def run_spades(read_dict):
    for key, value in read_dict.items():
        os.chdir(f"{key}") #navigate to folder containing the specific genome's read file(s)
        if value == "single":
            filename = f"{key}"+".fastq" #name of file containing reads for this genome
            spades_command = f"spades.py -s {filename} -o {key}_assembly" #building command to run spades and output in the same folder where the reads are contained
        elif value == "paired":
            file1 = f"{key}"+"_1.fastq" #names of files containing reads for this genome
            file2 = f"{key}"+"_2.fastq"
            spades_command = f"spades.py -1 {file1} -2 {file2} -o {key}_assembly" #building command to run spades
        os.system(spades_command) #execute command
        os.chdir("..") #navigate back to "Fastqs" folder

print(run_spades(read_types))
print("SPAdes complete")
