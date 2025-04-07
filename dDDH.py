import os
import csv

os.chdir("Ecoli-Project1")

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
                with open(f"{value[i]}_{value[i+j]}_blast.tsv") as file: #access blast results
                    rd = csv.reader(file, delimiter="\t")
                    lengthlist = []
                    alignlist = []
                    nidentlist = []
                    for row in rd:
                        nident = int(row[0])
                        nidentlist.append(nident)
                        length = int(row[1])
                        lengthlist.append(length)
                        qcovus = (int(row[2]))/100
                        adjust_align = nident * qcovus #adjusted alignment length
                        alignlist.append(adjust_align)
                    # print(alignlist)
                    # print(lengthlist)
                regdDDHlist = []
                adjdDDHlist = []
                for x in range(len(lengthlist)): #calculate regular and adjusted dDDH values for each contig
                    regdDDHlist.append(nidentlist[x] / lengthlist[x])
                    adjdDDHlist.append(alignlist[x] / lengthlist[x])
                # calculating average dDDH across all contigs for both regular and adjusted dDDH
                regdDDH = 0
                for item in regdDDHlist:
                    regdDDH += item
                regdDDH = (regdDDH / len(regdDDHlist))
                adjdDDH = 0
                for item in adjdDDHlist:
                    adjdDDH += item
                adjdDDH = (adjdDDH / len(adjdDDHlist))   
                print(regdDDH)
                print(adjdDDH)             
                # print(regdDDHlist)
                # print(adjdDDHlist)
                # for y in range(len(regdDDHlist)):
                #     print(regdDDHlist[y], adjdDDHlist[y])
                # print()

            # blast genome pairs and output length, qcovus, and nident
            # dDDH = (nident * qcovus) / length
            # store dDDH in tsv with first column of genome vs genome
    return("deez")

print(run_blast(groups))