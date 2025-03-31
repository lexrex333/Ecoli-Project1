import os
import csv

#create a tsv file to hold output data that we can work with
f = open("output_data.tsv", "w")
f.close()

# outfmt length nident
# Alignment length, num of identical matches
# sum of all identities found in matches divided by the overall match length = dDDH

comparedict = {}

for key, value in comparedict.items():
    for genome in value:
        print(genome)
        # blast genome pairs and output length and nident
        # dDDH = nident / length
        # store dDDH in tsv with first column of genome vs genome