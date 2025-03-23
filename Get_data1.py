import csv
import subprocess
import os
os.system("mkdir Ecoli-Project1")
os.chdir("Ecoli-Project1")
#empty list to hold all the accession numbers
accessions = []

main_dir = "/home/project1/Ecoli-Project1" #main directory - change here 
os.makedirs(main_dir, exist_ok=True) #make the directory if not there
os.chdir(main_dir) #move into the directory

csv_strains = os.path.join(main_dir, "LTEE_strains.csv") #csv with strain info
sra_prefetch_output = os.path.join(main_dir, "Prefetch_files") #where the sra files are at
fastqs_faster_output = os.path.join(main_dir, "Fastqs") #where fasterq output is going to go

os.makedirs(sra_prefetch_output, exist_ok=True) #make directory if not there
os.makedirs(fastqs_faster_output, exist_ok=True) #same here - make directory if not there


# Step 1: opening csv file with all the strain info and moving the SRR numbers to a list
with open(csv_strains, mode = 'r') as f: #opening and indicating reading the csv file
    csv_dict = csv.DictReader(f) #creating a dictionary for each row where the headers are the keys
    for row in csv_dict: #for each line in the csv dictionary
        accession = row["Accession"].strip() #getting the accession numbers 
        accessions.append(accession) #and adding it to assessions list
        #print(row)
    
#printing to see that the accession numbers were extracted correctly 
print(accessions)


# Step 2: get fastqs from NCBI
for accession in accessions: #for each accession in the accessions list
    print(type(accession))
    #prefetch gets the SRA number from NCBI
    prefetch = ["prefetch", "-t", "2", accession, "-O", sra_prefetch_output] #outputting into fastqs_NCBI directory
    print(prefetch)
    subprocess.run(prefetch, check = True) #running the prefetch command

    #go to path with sra file that was downloaded by prefetch
    sra_pathway = os.path.join(sra_prefetch_output, accession, f"{accession}.sra")
    #make directories for fasterq output 
    output_dir = os.path.join(fastqs_faster_output, accession)
    os.makedirs(output_dir, exist_ok=True)

    #running fasterq to make the 2 paired reads for each sra number
    fasterq = ["fasterq-dump", "-t", "2", sra_pathway, "-O", output_dir]
    #running the fasterq command
    subprocess.run(fasterq, check = True)

print("All fastqs have been downloaded.")


  