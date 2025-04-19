import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

#won't need this in combined script
main_dir = "/home/2025/aavalos4/Ecoli-Project1"

#getting files and making dataframes
#ani/ other number results 
num_df = pd.read_csv(os.path.join(main_dir, "FastANI_Results_ALL.csv"), header = None)

#where all the info is at
all_info = pd.read_csv(os.path.join(main_dir, "LTEE_strains.csv"))

#putting columns to the ani dataframe
num_df.columns = ["Query", "Reference", "ANI", "Gen1","Gen2", "Population", "Label"]

#want only the comparisons to the ancestor rows
num_df = num_df[num_df["Label"] == "To Ancestor"]

#keeping the columns we need to compare
num_df = num_df[["Query","ANI"]]

#merging the data
merge_df = num_df.merge(all_info, left_on = "Query", right_on= "Accession")

#keeping only the columns we need
merge_df =merge_df[["Accession", "Generation", "ANI"]]

#debugging - changing the generation to make sure its in number format
merge_df["Generation"]= pd.to_numeric(merge_df["Generation"])
merge_df["ANI"]= pd.to_numeric(merge_df["ANI"])

#sorting by the generations
merge_df = merge_df.sort_values("Generation")

#creating data for heatmap format - one row matrix for ani values bc one pop
heatmap_data = pd.DataFrame([merge_df["ANI"].values],columns=merge_df["Generation"].values) 

#plotting
plt.figure(figsize=(15,2)) #figure size
sns.heatmap(heatmap_data, cmap="viridis", cbar_kws={"label":"ANI to Ancestor"})
plt.xticks(rotation=90) #x-axis
plt.yticks([],[]) #hiding y-axis
plt.xlabel("Generation") #x-axis
plt.title("ANI to Ancestor by Generation") #title of plot
plt.tight_layout() #making sure the layout is nice
plt.savefig(os.path.join(main_dir,"ANI_Heatmap.png"), dpi=300, bbox_inches="tight")#saving to output file
plt.show()

#for group code
#1. need to manually download seaborn 
#2. need to make the code doable for populations and not just one population
