import csv
import subprocess
import os

#empty list to hold all the accession numbers
accessions = []

# Step 1: opening csv file with all the strain info and moving the SRR numbers to a list
with open('/home/project1/LTEE_strains.csv', mode = 'r') as f: #opening and indicating reading the csv file
    csv_dict = csv.DictReader(f) #creating a dictionary for each row where the headers are the keys
    for row in csv_dict: #for each line in the csv dictionary
        accessions.append(row["Accession"].strip()) #getting the accession numbers and adding it to assessions list
        #print(row)
    
#printing to see that the accession numbers were extracted correctly 
print(accessions)


# Step 2: get fastqs from NCBI
for accession in accessions: #for each accession in the accessions list
    #prefetch gets the SRA number from NCBI
    prefetch = ["prefetch", "-t", 2, accession, "-O", "home/project1/fastqs_NCBI" ] #outputting into fastqs_NCBI directory
    subprocess.run(prefetch, check = True) #running the prefetch command
    
os.system("cd fastqs_NCBI") #going into the directory

# Step 3: changing sra numbers from NCBI into paired reads fastqs
for num in accessions:
    
    #running fasterq to make the 2 paired reads for each sra number
    fasterq = ["fasterq-dump", "-t", 2, num ]


