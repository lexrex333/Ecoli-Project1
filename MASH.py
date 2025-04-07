import os
import csv
accessions = []
groups = {
    "Ara+1": []
}

with open("sample_data.csv", mode = 'r') as f: #opening and reading the csv file
	csv_dict = csv.DictReader(f) #creating a dictionary for each row where the headers are the keys
	for row in csv_dict: #for each line in the csv dictionary
		accession = row["Accession"].strip() #getting the accession numbers and stripping extra whitespace
		accessions.append(accession) #and adding it to accessions list
		group = row["\ufeffPopulation"].strip()
		if group in groups:
			groups[group].append(accession)
		for group in groups:
			groups[group].append("SRR22764941") #add ancestor accession to the end of each population group
for num in accessions:
	if num != "SRR22764941":
		fastq_file = 
def MASH(one,two):
	mash_command = ["mash dist ", one, two]
	print(mash_command)
	subprocess.run(mash_command)
