import os
import subprocess

# Paths and parameters; Edit them as needed. No input upon running the script, just make sure to set the paths correctly before execution. @GmoreiraVet
input_folder = "/home/viroicbas2023/Downloads/Bluetongue/JTF5YV_fastq"
output_folder = "/home/viroicbas2023/Downloads/Bluetongue/Reports"
kraken_db = "/home/viroicbas2023/Documents/Gmoreira/krakenDB/PlusPfP_GrandalhonaBro"

# Toggle exporting of unclassified reads. Potentially useful for uncharacterized sequences; Mostly useful for trash reads. True or Flase. @GmoreiraVet
EXPORT_UNCLASSIFIED = True

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Kraken2 parameters.
# Default Parameters: --threads 30 (for Watson & Crick only, check your system for available threads), --minimum-base-quality 8, --confidence 0, --memory-mapping
threads = "30"
min_base_quality = "8"
confidence = "0"
memory_mapping = "--memory-mapping"

# Iterate over all FASTQ files in the input folder
#os.listdir() lists all files in a folder.
#File name is checked for FASTQ extensions (.fq .fastq). Only these will be affected.
for file_name in os.listdir(input_folder):
    if file_name.endswith(".fastq") or file_name.endswith(".fq"):
        #Input_path function constructs the full path by joining folder path with file name.
        input_path = os.path.join(input_folder, file_name)
        #Extracts the base name of the file (without extension) to use in output. [0] is used to get the first part of the split, which is the base name without extension.
        sample_name = os.path.splitext(file_name)[0]
        # Output paths for Kraken2 results and reports, constructed using the sample name.
        output_path = os.path.join(output_folder, f"{sample_name}_kraken.out")
        report_path = os.path.join(output_folder, f"{sample_name}_kraken.report")

        # Base Kraken2 command. The command is built as a list of arguments, which is the recommended way to use subprocess.run() to avoid issues with spaces and special characters in file names. Each parameter is added to the command list, and the input file is appended at the end. If EXPORT_UNCLASSIFIED is True, additional arguments are added to export unclassified reads to a separate FASTQ file. Finally, the command is executed using subprocess.run(), which runs Kraken2 with the specified parameters for each FASTQ file in the input folder.
        cmd = [
            "kraken2",
            "--db", kraken_db,
            "--threads", threads,
            "--output", output_path,
            "--report", report_path,
            "--minimum-base-quality", min_base_quality,
            memory_mapping,
            "--confidence", confidence,
        ]

        # Optional export of unclassified reads. If EXPORT_UNCLASSIFIED is set to True, the script constructs an additional output path for unclassified reads and appends the necessary arguments to the Kraken2 command to export these reads to a separate FASTQ file. This allows users to retain unclassified sequences for further analysis or quality control, which can be particularly useful for identifying novel or poorly characterized organisms in metagenomic samples.
        if EXPORT_UNCLASSIFIED:
            unclassified_out = os.path.join(output_folder, f"{sample_name}_unclassified.fastq")
            cmd.extend(["--unclassified-out", unclassified_out])

        # Input file. The input FASTQ file is appended to the command list, which is the last argument for Kraken2. This ensures that Kraken2 processes the correct file for each iteration of the loop, allowing for efficient batch processing of multiple samples in the specified input folder.
        cmd.append(input_path)

        # Run Kraken2 command. The script prints a message indicating which file is being processed, executes the Kraken2 command using subprocess.run(), and then prints a message upon completion of each file, including the paths to the output and report files. This provides feedback to the user about the progress of the analysis and allows for easy tracking of results.
        print(f"Running Kraken2 on {file_name}...")
        subprocess.run(cmd)
        print(f"Finished {file_name}. Output: {output_path}, Report: {report_path}\n")

print("All samples processed.")
