import os
import csv
import subprocess
from itertools import combinations

accessions = []
groups = {
    "Ara+1": []
}

base_path = "/home/2025/jfloros/Comp_Bio/FinalProject/sample_assemblies"
results_file = "mash_results.txt"

# Clear the file at the beginning (optional, just to start fresh each run)
with open(results_file, "w") as f:
    f.write("MASH Distance Results\n")

# Read CSV
with open("sample_data.csv", mode='r') as f:
    csv_dict = csv.DictReader(f)
    for row in csv_dict:
        accession = row["Accession"].strip()
        accessions.append(accession)
        group = row["\ufeffPopulation"].strip()
        if group in groups:
            groups[group].append(accession)

# Add ancestor accession
ancestor = "SRR22764941"
for group in groups:
    groups[group].append(ancestor)

# Get FASTQ file path
def get_contigs_path(accession):
    return os.path.join(base_path, f"{accession}_assembly", "contigs.fasta")

# Create a sketch for an accession, return the .msh file path
def create_sketch(contigs_path):
    sketch_path = contigs_path + ".msh"
    if os.path.exists(sketch_path):  # Force fresh sketch with updated paths
        os.remove(sketch_path)
    #if not os.path.exists(sketch_path):  # Only sketch if it doesn't exist
    sketch_command = ["./mash", "sketch", "-o", sketch_path, contigs_path]
    subprocess.run(sketch_command, check=True)
    return sketch_path

# Run MASH and write results
def MASH(one_accession, two_accession):
    one_contigs = get_contigs_path(one_accession)
    two_contigs = get_contigs_path(two_accession)

    sketch_one = create_sketch(one_contigs)
    sketch_two = create_sketch(two_contigs)

    mash_command = ["./mash", "dist", sketch_one, sketch_two]
    result = subprocess.run(mash_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    with open(results_file, "a") as f:
        f.write(f"{one_accession} vs {two_accession}:\n")
        f.write(result.stdout)
        if result.stderr:
            f.write("Errors:\n" + result.stderr)
        f.write("\n")

# Loop through pairwise comparisons
for group_name, accession_list in groups.items():
    for one, two in combinations(accession_list, 2):
        MASH(one, two)
