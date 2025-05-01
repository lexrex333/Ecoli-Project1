import csv
import subprocess
import os


main_dir = os.getcwd() #putting variable for the directory you are in

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

