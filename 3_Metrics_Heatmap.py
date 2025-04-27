import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os


#will change when use 
main_dir = "/home/2025/aavalos4/Ecoli-Project1"

#Step 1. Read the tsv file of the 3 metrics 
ani = pd.read_csv(os.path.join(main_dir, "full_output_data.tsv"), sep="\t") #tab separated

#renaming the columns to match what I got going on
ani = ani.rename(columns={
    "sample_1": "Query", #the first strain
    "sample_2": "Reference", #the strain it is being compared to - ancestor
    "population": "Group" #the population it is comign from 
})

#Step 2. CLeaning everything
#Keeping only the comments I need and want in the heatmap
summary = ani[["Group", "ANI", "MASH_sim", "dDDH"]].copy()

#grouping by population 
pop_summary = summary.groupby("Group").mean()

#multiply MASH and dDDH to match the ani numbers 
pop_summary["MASH_sim"] = pop_summary["MASH_sim"] * 100 
pop_summary["dDDH"] = pop_summary["dDDH"] * 100

#Putting the populations in order - alphabetically
pop_summary = pop_summary.sort_index()

#Step 3. PLotting the heatmap
plt.figure(figsize=(8, 10)) #figure - size  

#plotting the heatmap 
sns.heatmap( #using seaborn 
    pop_summary,
    cmap="viridis", #the colors for the heatmap
    annot=True, #showing the number inside each square 
    fmt=".3f", #3 decimal point for the numbers 
    cbar_kws={"label": "Similarity Percentage"} #label for colorbar
)

#titles and axis labels 
plt.title("ANI, MASH Similarity, and dDDH Compared to Ancestor", fontsize=16) 
plt.ylabel("Population", fontsize=14) #y-axis label
plt.xlabel("Metric", fontsize=14) #x-axis label

#making the layout clean
plt.tight_layout()

#saving the figure 
out_path = os.path.join(main_dir, "3_Metric_Heatmap.png")
plt.savefig(out_path, dpi=300) #saving the figure 
plt.close() #closing
