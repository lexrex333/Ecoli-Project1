def run_trimmomatic():
    fastqs_dir= os.path.join(main_dir, "sample_reads") #set directory containing raw fastq files  
    trim_output= os.path.join(main_dir, "Trimmomatic_Results") #set directory to store the trimmed output files
   
    os.makedirs(trim_output, exist_ok=True) #make sure the output directory exsists                                    
    # Define the paths to the trimmomatic jar and adapter file 
    trimmomatic_jar = os.path.join(main_dir, "Trimmomatic-0.39", "trimmomatic-0.39.jar")
    adapter_file = os.path.join(main_dir, "Trimmomatic-0.39", "adapters", "TruSeq3-SE.fa")
    
    for root, dirs, files in os.walk(fastqs_dir): #go through the fastq directory and subdirectories to find file 
     for file in files:  #loop through all files ends with .fastq
        if file.endswith(".fastq"):
            fastq_path = os.path.join(root, file)  #create full path to the file
            #construct the path for the trimmed output file
            trimmed_fastq = os.path.join(trim_output, f"{file}")

            #construct the Trimmomatic command for trimming the fastq files 
            #run trimmomatic, single end reads, specify phred quality score, input Fastq file path and output trimmed file path
            #remove leading and trailing low quality bases, trimming if average quality <15, only retain reads that are 36 bases or longer 
            trimmomatic_command = [
                "/home/project1/miniconda3/envs/trimmomatic_env/bin/trimmomatic" , "SE", "-phred33", fastq_path, trimmed_fastq,
                 f"ILLUMINACLIP:{adapter_file}:2:30:10", "LEADING:3", "TRAILING:3",
                "SLIDINGWINDOW:4:15", "MINLEN:36"
            ] 
            print("Running:", " ".join(trimmomatic_command))
 #run the Trimmomatic command
            subprocess.run(trimmomatic_command, check=True)

     #print a message indicating that Trimmomatic processing is complete for all files
    print("Trimmomatic is done.")
    return trim_output

trim_output = run_trimmomatic() #running trimmomatic
 