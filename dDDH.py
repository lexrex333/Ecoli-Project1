import os
import csv

# accessions = []
groups = {
    "Ara+1": []
}

with open("sample_data_copy.csv", mode = 'r') as f: #opening and reading the csv file
    csv_dict = csv.DictReader(f) #creating a dictionary for each row where the headers are the keys
    for row in csv_dict: #for each line in the csv dictionary
        accession = row["Accession"].strip() #getting the accession numbers and stripping extra whitespace
        # accessions.append(accession) #and adding it to accessions list
        group = row["\ufeffPopulation"].strip()
        if group in groups:
            groups[group].append(accession)
        # print(row)
    for group in groups:
        groups[group].append("SRR22764941") #add ancestor accession to the end of each population group
# print(accessions)
print(groups)

#create a tsv file to hold output data that we can work with
# f = open("output_data.tsv", "w")
# # create a header row probably
# f.close()

# outfmt length qcovus nident
# Alignment length, query coverage without overlaps, num of identical matches
# sum of all identities found in matches divided by the overall match length = dDDH

os.chdir("sample_reads")

def run_blast(comparedict):
    for key, value in comparedict.items():
        for i in range(len(value)):
            if i+1 == len(value):
                continue
            for j in range(1, len(value)):
                if i+j >= len(value):
                    continue
                print(value[i], value[i+j])
                blast_command = f"blastn -query {value[i]}/{value[i]}_assembly/contigs.fasta -subject {value[i+j]}/{value[i+j]}_assembly/contigs.fasta -out {value[i]}_{value[i+j]}_blast.tsv -outfmt '6 nident length qcovus'"
                # os.system(blast_command)
                with open(f"{value[i]}_{value[i+j]}_blast.tsv") as file:
                    rd = csv.reader(file, delimiter="\t")
                    nident = []
                    adjlength = []
                    for row in rd:
                        nident.append(int(row[0]))
                        length = int(row[1])
                        qcovus = (int(row[2]))/100
                        adjlength.append(length * qcovus)
                print(nident)
                print(adjlength)
                dDDHlist = []
                for x in range(len(nident)):
                    dDDHlist.append(nident[x] / adjlength[x])
                print(dDDHlist)


            # blast genome pairs and output length, qcovus, and nident
            # dDDH = nident / (length * qcovus)
            # store dDDH in tsv with first column of genome vs genome
    return("deez")

print(run_blast(groups))