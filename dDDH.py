import os
import csv

os.chdir("Ecoli-Project1")

groups = {
    "Ara+1": []
}

with open("sample_data_copy.csv", mode = 'r') as f: #opening and reading the csv file
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
    f.close()
    for key, value in comparedict.items():
        for i in range(len(value)):
            if i+1 == len(value):
                continue
            for j in range(1, len(value)):
                if i+j >= len(value):
                    continue
                print(value[i], value[i+j])
                blast_command = f"blastn -query sample_assemblies/{value[i]}_assembly/contigs.fasta -subject sample_assemblies/{value[i+j]}_assembly/contigs.fasta -max_hsps 1 -out blast_output/{value[i]}_{value[i+j]}_blast.tsv -outfmt '6 nident length'"
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
                f = open("blast_output_data.tsv", "w")
                f.write(f"{value[i]}\t{value[i+j]}\t{dDDH}\n")
                f.close()
    return("deez")

print(run_blast(groups))