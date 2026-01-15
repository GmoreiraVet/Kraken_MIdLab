import os
import subprocess

# Paths
input_folder = "/home/viroicbas2023/Documents/Gmoreira/ChildSJoao/MXN8CN_fastq"
output_folder = "/home/viroicbas2023/Documents/Gmoreira/ChildSJoao/KrakenHugeConf01"
kraken_db = "/home/viroicbas2023/Documents/Gmoreira/krakenDB/PlusPfP_GrandalhonaBro"

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Kraken2 parameters
threads = "30"
min_base_quality = "8"
confidence = "0.1"
memory_mapping = "--memory-mapping"

# Iterate over all FASTQ files in the input folder
for file_name in os.listdir(input_folder):
    if file_name.endswith(".fastq") or file_name.endswith(".fq"):
        input_path = os.path.join(input_folder, file_name)
        sample_name = os.path.splitext(file_name)[0]
        output_path = os.path.join(output_folder, f"{sample_name}_kraken.out")
        report_path = os.path.join(output_folder, f"{sample_name}_kraken.report")

        # Build Kraken2 command
        cmd = [
            "kraken2",
            "--db", kraken_db,
            "--threads", threads,
            "--output", output_path,
            "--report", report_path,
            "--minimum-base-quality", min_base_quality,
            memory_mapping,
            "--confidence", confidence,
            input_path
        ]

        print(f"Running Kraken2 on {file_name}...")
        subprocess.run(cmd)
        print(f"Finished {file_name}. Output: {output_path}, Report: {report_path}\n")

print("All samples processed.")
