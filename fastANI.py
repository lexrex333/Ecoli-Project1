import os
import subprocess
import pandas as pd
from itertools import combinations



main_dir = '/home/2025/aavalos4/Ecoli-Project1' #remove this when adding to pipeline



#this function is going to run fastani within each pop group and compares each sample to the ancestor strain
def fastani_run(main_dir, sample_reads, sample_csv, fastani_output="FastANI_Results", ancestor_name = "REL606"):
    
    
    #THIS SECTION Will be deleted --- can grab this from previous code place where the sample reads are at
    sample_reads = os.path.join(main_dir, 'sample_reads') 
    #place to csv with data
    sample_csv = os.path.join(main_dir, "sample_data.csv")
    #--------------------

    #place to save fastani output
    fastani_output = os.path.join(main_dir, 'FastANI_Results')
    #make folder if not already made
    os.makedirs(fastani_output, exist_ok=True)

    #loading sample data to look off of
    df = pd.read_csv(sample_csv)

    sample_fastas = {} #place to put paths of each spades fasta file

    #dictionary to hold the population groups
    pop_groups = {}

    #going through each row of the csv to group the fastas based on their groups
    for _, row in df.iterrows():
        sample = row["Accession"].strip() #getting accession row
        population = row["Population"].strip() #getting population row
        
        #path to spades fastas
        fasta_path = os.path.join(main_dir,"sample_reads", f"{sample}_assembly", "contigs.fasta") #fasta path to spades output file
        
        if os.path.isfile(fasta_path): #if there actually is a file in that path
            sample_fastas[sample] = fasta_path #add the folder name and path to the fasta file to the fastas dictionary 

            #adding the sample to the pop group that it belongs to
            if population not in pop_groups: #if the population is already not added to the groups dictionary
                pop_groups[population] = [] #add the population as a key
            pop_groups[population].append(sample) #add the sample to the population 
        
        #debugging 
        else:
            print(f"No fasta for {sample} at {fasta_path}")

    #Step ONE: comparing within each population
    for population, samples in pop_groups.items():
        #making a folder for just this population
        this_pop = os.path.join(fastani_output, f"{population}_within")
        os.makedirs(this_pop, exist_ok=True) #add if not there already

        #getting all the comparisons bw 2 samples between this individual group
        #combinations(samples,2) is going to do all unique pairs without repeating them or comparing them to themselves
        for sample1, sample2 in combinations(samples, 2):
            fasta1 = sample_fastas[sample1]
            fasta2 = sample_fastas[sample2]

            #naming the output file for the pairwise comparison 
            output_fastANI = os.path.join(this_pop, f"{sample1}_vs_{sample2}.txt") #each file is going to be named "fasta_vs_fasta.txt" to show what was compared

            #fastani command
            fastani = ["fastANI", "--query", fasta1, "--reference", fasta2, "--output", output_fastANI]
            #running the command
            subprocess.run(fastani, check=True)

    #STEP TWO: Comparing ALL the samples to the ancestor
    ancestor = sample_fastas.get(ancestor_name) #might have to change this later when added to pipeline

    if not ancestor: #if ancestor name not found in the sample list
        return
    
    #making folder for the comparisons to the ancestor outputs to go
    ancestor_output = os.path.join(fastani_output, f"{ancestor_name}_comparisons")
    os.makedirs(ancestor_output, exist_ok=True) #make if not there already

    #comparing each fasta to the ancestor except the ancestor fasta
    for sample, fasta in sample_fastas.items(): #going through sample_fastas dictionary
        if sample == ancestor_name: #if the sample equals to the ancestor
            continue #ignore and keep going

        #making path to put the output txt files
        ancestor_compare_output = os.path.join(ancestor_output, f"{sample}_vs_{ancestor_name}.txt") #naming each output file based on names for clarity
        #the fastani command
        fastani_command = ["fastANI", "-q", fasta, "-r", ancestor, "-o", ancestor_compare_output]
        #running the fastani command
        subprocess.run(fastani_command, check = True)
        
    print("FastANI has been completed.")

fastani_run(main_dir, sample_csv, fastani_output = "FastANI_Results", ancestor_name = "insert sample name")

#need to make sure the variables match up
#also might need to write to one file rather than multiple txt files  
#maybe make function for the population grouping rather than just have it in this function?