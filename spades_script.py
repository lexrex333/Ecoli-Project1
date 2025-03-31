import os
import csv

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