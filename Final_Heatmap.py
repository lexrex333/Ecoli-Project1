import pandas as pd  
import seaborn as sns   
import matplotlib.pyplot as plt 
import os    

#wont need this later in final code
main_dir = "/home/2025/aavalos4/Ecoli-Project1"

#putting the ancestor accession so it can identify it in the heatmap
ancestor_accession = "SRR22764941" #ancestor accession
ancestor_label = f"{ancestor_accession} (ancestor)" #when plotting, will give accession number of ancestor with ancestor next to it

#Step 1: getting the data for strains
#getting the fastani or other metric results
ani = pd.read_csv(os.path.join(main_dir, "FastANI_Results_ALL.csv"), header=None) #will have to change this for another metric
#columns found in the fastani result columns 
ani.columns = ["Query", "Reference", "ANI", "Matching_Fragments", "Total_Fragments", "Group", "Comparison_Method"] #might have to change this for other metric csv to match

#getting the strain info from the original csv with all strain info
strain_info = pd.read_csv(os.path.join(main_dir, "LTEE_strains.csv"))


#Step 2: cleaning up so i dont get an error
#for each of the columns 
for col in ["Query", "Reference", "ANI", "Group"]: #that are identified here
    ani[col] = ani[col].astype(str).str.strip() #changing to string and taking out whitespace

#making sure the ani values are numbers to avoid error -debugging
ani["ANI"] = pd.to_numeric(ani["ANI"], errors="coerce") #if cant be turned into a number, 1.red flag, 2. it will put NAN and we can see somehting is wrong

#also doing the same thing for strain info columns 
for col in ["Population", "Generation", "Accession"]:
    strain_info[col] = strain_info[col].astype(str).str.strip() #string and cleaning whitespace

#cnanging the generation to an integer, and making a new label with the SRR num and generation
strain_info["Generation"] = pd.to_numeric(strain_info["Generation"], errors="coerce").astype("Int64") #will support NAN wothout quitting out
strain_info["Label"] = strain_info["Accession"] + " (Gen " + strain_info["Generation"].astype(str) + ")" #new name/label here 


#Step 3: want to add generation and label info to ani dataframe
#so we know what info goes to each strain -info we want to use
#so want to create dictionaries that connect each SRR num to its matching data

#so this makes a dictionary where SRR#-> SRR# (Gen #)
label_map = dict(zip(strain_info["Accession"], strain_info["Label"]))

#this makes a dictionary where ex: SRR#-> 500 
gen_map = dict(zip(strain_info["Accession"], strain_info["Generation"]))

#this makes a dictionary where SRR#-> Ara+1 (population)
pop_map = dict(zip(strain_info["Accession"], strain_info["Population"]))

#now we want to use these dictionaries to put stuff on the ani dataframe 
#for each query strain, add its real label we want on the heatmap
ani["Query_Label"] = ani["Query"].map(label_map)

#for each reference strain, do the exact same thing we did for the query strain
ani["Reference_Label"] = ani["Reference"].map(label_map)

#adding the generation number for all the strains 
ani["Query_Gen"] = ani["Query"].map(gen_map) #query strains
ani["Ref_Gen"] = ani["Reference"].map(gen_map) #reference strains

#adding the pop info for all the strains
ani["Query_Pop"] = ani["Query"].map(pop_map) #query strains
ani["Ref_Pop"] = ani["Reference"].map(pop_map) #reference strains


#Step 4: now making the heatmaps for each pop
#go through each of the pops separately
for population in ani["Query_Pop"].dropna().unique(): #removes missing values, returning only unique vlaues
    
    #making a row to represent the ancestor
    ancestor_row = pd.DataFrame({
        "Accession": [ancestor_accession], #accession number
        "Generation": [0], #because it started at the ancestor
        "Population": [population], #will match current pop, but should be the same for all pops
        "Label": [ancestor_label] #putting in the label for the heatmap we want
    })

    #adding this just in case its not there - debugging 
    if not strain_info[strain_info["Accession"] == ancestor_accession].any().any():
        strain_info = pd.concat([strain_info, ancestor_row], ignore_index=True) #adding if not there

    #only keeping the ani rows that both the query and reference strains belong to the pop
    pop_df = ani[(ani["Query_Pop"] == population) & (ani["Ref_Pop"] == population)].copy()

    #creaitng a grid [pivot table]
    ani_table = pop_df.pivot_table(
        index="Query_Label", #rows = query strains
        columns="Reference_Label", #columns = reference strains
        values="ANI", #values = ani scores 
        aggfunc="mean" #if there is more than one value, will take the average, but there shouldn't be bc each SRR is only in there once 
    )

    #FastAni only compares in one direction so SRR -> other SRR, and not other SRR -> SRR
    # so we need to mirror them to fill in the heatmap
    #.T = transpose -- helps flip it to fill in the other side
    ani_table = ani_table.combine_first(ani_table.T)

    #Step 5: sorting the heatmap by generation so its easier to read
    #taking out the generation num from the labels
    def extract_gen(label):
        if "ancestor" in label: #if the ancestor is in the label
            return 0 #we don't want it bc its already at 0
        elif "Gen " in label: #if there is gen in there 
            return int(label.split("Gen ")[1].replace(")", "")) #then most likely SRR# strain from pops and want to grab the gen number 
        else: #anything else
            return 99999  #so we know something went wrong - debugging

    #now sorting both the rows and columns by the generation using function above - sorting from oldest to most recent
    ani_table = ani_table.reindex(sorted(ani_table.index, key=extract_gen)) #sorting rows
    ani_table = ani_table[sorted(ani_table.columns, key=extract_gen)] #sorting columns

    #Step 6: making and saving heatmap

    plt.figure(figsize=(12, 10)) #size of plot

    #making the heatmap 
    sns.heatmap(
        ani_table, #using our dataframe
        cmap="viridis", #color for heatmap
        linewidths=0.5, #lines between all the boxes
        linecolor="gray", #color of lines
        vmin=85, vmax=100, #the range of the ani values on color scale
        cbar_kws={"label": "ANI"} #labeling on the color bar -- change for other metric
    )

    #adding title and axes labels
    plt.title(f"Pairwise ANI Heatmap - Population {population}")
    plt.xlabel("Strain (Accession + Generation)")
    plt.ylabel("Strain (Accession + Generation)")

    #making it so we can read the labels better
    plt.xticks(rotation=90) #x-axis
    plt.yticks(rotation=0) #y-axis
    plt.gca().invert_yaxis() #to make the y-axis go the opposite way so increasing upward instead of decreasing upward
    plt.tight_layout() #finalizing layout

    #saving the heatmap to a png file with pop name - will make a new one for each pop
    out_path = os.path.join(main_dir, f"ANI_Heatmap_{population}.png") #place where to put png 
    plt.savefig(out_path, dpi=300) #saving the figure
    plt.close() 
