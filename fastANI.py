import os
import subprocess
import pandas as pd
from itertools import combinations



main_dir = '/home/2025/aavalos4/Ecoli-Project1' #remove this when adding to pipeline



#FastANI

#this function is going to run fastani within each pop group and compares each sample to the ancestor strain
def fastani_run(main_dir, sample_csv, fastani_output="FastANI_Results", ancestor_name = "REL606A.fasta"):

    #place to save fastani output
    fastani_output = os.path.join(main_dir, 'FastANI_Results')
    #make folder if not already made
    os.makedirs(fastani_output, exist_ok=True)

    #loading sample data to look off of
    df = pd.read_csv(sample_csv)

    sample_fastas = {} #place to put paths of each spades fasta file

    #dictionary to hold the population groups
    pop_groups = {}

    #place to put all FastANI results before putting it to csv file
    results = [] #empty list

    #going through each row of the csv to group the fastas based on their groups
    for _, row in df.iterrows(): #_ ignores the index and just gets the rows
        sample = row["Accession"].strip() #getting the sample name
        population = row["Population"].strip() #getting population row
        
        #path to spades fastas
        fasta_path = os.path.join(main_dir,"sample_assemblies", f"{sample}_assembly", "contigs.fasta") #fasta path to spades output file
        
        if os.path.isfile(fasta_path): #if there actually is a file in that path
            sample_fastas[sample] = fasta_path #add the folder name and path to the fasta file to the fastas dictionary 

            #adding the sample to the pop group that it belongs to
            if population not in pop_groups: #if the population is already not added to the groups dictionary
                pop_groups[population] = [] #add the population as a key
            pop_groups[population].append(sample) #add the sample to the population 
        
        #debugging 
        else:
            print(f"No fasta for {sample} at {fasta_path}")

    #debugging
    print(list(sample_fastas.keys()))
    print(pop_groups.keys())

    #Step ONE: comparing within each population
    for population, samples in pop_groups.items():

        #getting all the comparisons bw 2 samples between this individual group
        #combinations(samples,2) is going to do all unique pairs without repeating them or comparing them to themselves
        for sample1, sample2 in combinations(samples, 2):
            fasta1 = sample_fastas[sample1]
            fasta2 = sample_fastas[sample2]

            #fastani command -- will need to change the pathway by manually getting fastANI
            fastani = [main_dir, "fastANI", "-q", fasta1, "-r", fasta2, "--stdout"] #storing output to terminal instead to write to csv
            #running the command
            result = subprocess.run(fastani, capture_output=True, text=True, check=True)
            #print(result)
            
            #going through the output and adding to the results list
            for line in result.stdout.strip().splitlines():
                _, _, ani, fragments, total = line.strip().split("\t") #strip whitespace and split by tab
                results.append({
                    "Query": sample1,
                    "Reference": sample2,
                    "ANI": ani,
                    "Matching_Fragments": fragments,
                    "Total_Fragments": total,
                    "Group": population,
                    "Comparison_Method": "Within Population"
                })

    #STEP TWO: Comparing ALL the samples to the ancestor
    ancestor = sample_fastas.get(ancestor_name) #might have to change this later when added to pipeline

    if not ancestor: #if ancestor name not found in the sample list
        #debugging
        print(f"Ancestor {ancestor_name} not found.")
        return
    
    #comparing each fasta to the ancestor except the ancestor fasta
    for sample, fasta in sample_fastas.items(): #going through sample_fastas dictionary
        if sample == ancestor_name: #if the sample equals to the ancestor
            continue #ignore and keep going

        #the fastani command
        fastani_command = [main_dir, "fastANI", "-q", fasta, "-r", ancestor, "--stdout"]
        #running the fastani command
        result2 = subprocess.run(fastani_command, capture_output=True, text = True, check = True) #cleaner way of getting results
        #print(result2)

        #this looks through the csv and gets the row where the accession is equal to the sample
        #and gets the value that is in the population column that matches that 
        #in short talk, it is just getting the pop group that the sample belongs to
        group = df[df["Accession"] == sample]["Population"].values[0]

        #going through the output and adding to the results list
        for line in result2.stdout.strip().splitlines():
            _, _, ani, fragments, total = line.strip().split("\t") #strip whitespace and split by tab
            results.append({"Query": sample,
                            "Reference": ancestor_name,
                            "ANI": ani,
                            "Matching_Fragments": fragments,
                            "Total_Fragments": total,
                            "Group": group,
                            "Comparison_Method": "To Ancestor"
                            })
    
    results_to_df = pd.DataFrame(results) #changing the results to a dataframe
    results_to_df.to_csv(os.path.join(fastani_output, "FastANI_Results_ALL.csv"), index = False)
    
    print("FastANI has been completed!!")

fastani_run(main_dir, sample_csv = os.path.join(main_dir, "sample_data.csv"), fastani_output = "FastANI_Results", 
            ancestor_name = "SRR22764941")