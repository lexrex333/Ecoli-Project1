import pandas as pd
import os

#change this for the path of your main directory
main_dir = "/home/project1/Ecoli-Project1"

#path to the csv file
csv = os.path.join(main_dir, "LTEE_strains.csv") #will have to change the name to match your csv file

#path to the directory with the fastqs
fastqs = os.path.join(main_dir,"Fastqs") 

#loading csv file into a pandas dataframe
df = pd.read_csv(csv, header = None)

#getting the accession numbers and read type (whether its single or paired)
accessions = df[[6, 7]].dropna() #picking them out and dropping rows or columns that are empty in csv file - debug
accessions.columns = ["accession", "read_type"] #the names of the columns

#cleaning it up so it doesn't give errors
accessions["accession"] = accessions["accession"].astype(str).str.strip().str.split('.').str[0] #cleaning up --debugging
accessions["read_type"] = accessions["read_type"].str.strip().str.lower() #cleaning up

#keeping track of the missing files in a list
missing_fastq = []

#going through all the accession numbers
for i, row in accessions.iterrows():
    accession = row["accession"] #going through the rows of accession numbers - grabbing accession number
    read_type = row["read_type"] #going through the rows of the read types- grabbing read type
    fastqs_path = os.path.join(fastqs, accession) #path where fastqs are stored

    #seeing if the fastq is actually there
    if not os.path.isdir(fastqs_path): #if folder is not there
        missing_fastq.append(f"{accession}: folder is missing") #add to missing file list

    #if the sample is single end, it will only look for one file
    if read_type == "single": #if the read type is single
        file = os.path.join(fastqs_path, f"{accession}.fastq") #path to the fastq
        if not os.path.isfile(file): #if file not there
            missing_fastq.append(f"{accession}: single-end file missing") #report it

    #if the sampled is paired reads, it will look for 2 files
    elif read_type == "paired": #if it is paired
        file1 = os.path.join(fastqs_path, f"{accession}_1.fastq") #path to 1 fastq
        file2 = os.path.join(fastqs_path, f"{accession}_2.fastq") #path to 2 fastq
        if not (os.path.isfile(file1)): #if file not there
            missing_fastq.append(f"{accession}:_1.fastq missing") #records it 
        if not (os.path.isfile(file2)): #if file not there
            missing_fastq.append(f"{accession}:_2.fastq missing") #records it

    else: #else the read type is not readable - not single or paired
        missing_fastq.append(f"{accession}: we don't know the read type: {read_type}")

#giving info on results
if missing_fastq:
    print(f"Missing these fastq files:")
    for m in missing_fastq: #printing the fastq if missing
        print(m)
else:
    print(f"We matched all fastq files!! :)")