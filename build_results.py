import csv

f = open("full_output_data.tsv", "w") #create an output file
f.close()

#make lists to hold all the important imformation from these files
sample1list = []
sample2list = []
ANIlist = []
MASHlist = []
dDDHlist = []
poplist = []

#read dDDH file and pull information
with open("dddh_output.tsv") as file: #access data from tsv file
    rd = csv.reader(file, delimiter="\t")
    next(file) #skip first line
    for row in rd:
        sample1list.append(f"{row[0]}")
        sample2list.append(f"{row[1]}")
        dDDHlist.append(f"{row[2]}")
        poplist.append(f"{row[3]}")

#read MASH file and pull information
with open("mash_results.tsv") as file:
    rd = csv.reader(file, delimiter="\t")
    next(file) #skip first line
    for row in rd:
        MASHlist.append(f"{row[2]}")

#read ANI file and pull information
with open("FastANI_Results/FastANI_Results_ALL.csv") as file:
    rd = csv.reader(file, delimiter=",")
    next(file) #skip first line
    for row in rd:
        ANIlist.append(f"{row[2]}")

#write info in new output file
f = open("full_output_data.tsv", "a")
f.write("sample_1\tsample_2\tANI\tMASH\tdDDH\tpopulation\n") #create header row
for i in range(len(ANIlist)): #put info for each pair in each row
    f.write(f"{sample1list[i]}\t{sample2list[i]}\t{ANIlist[i]}\t{MASHlist[i]}\t{dDDHlist[i]}\t{poplist[i]}\n")
