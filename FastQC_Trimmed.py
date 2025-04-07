#another function to run fastqc because it didn't want to work for trimmed files :(
def fastqc_trimmed(folder): 
    fastq_files = get_fastq(folder) #calling function to make a fastqs list

    #creating the output folder
    fastqc_output = os.path.join(folder, "FastQC_Results") #where the results will go -- within the trimmomatic results
    os.makedirs(fastqc_output, exist_ok=True) #making directory if not already there

    for fastq in fastq_files: #for each of the fastqs within the fastq files list
        fastq_path = os.path.join(folder, fastq) #path of where the fastqs are at 
        fastqc_command = ["fastqc", fastq_path, "-o", fastqc_output] #fastqc command
        subprocess.run(fastqc_command, check = True) #run fastqc command 
