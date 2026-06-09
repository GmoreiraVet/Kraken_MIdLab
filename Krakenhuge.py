import os
import subprocess

# Paths and parameters; Edit them as needed. No input upon running the script, just make sure to set the paths correctly before execution. @GmoreiraVet
input_folder = "/home/viroicbas2023/Documents/Gmoreira/Carraças_Metagenomica_Gui/Reads_DeRNAfied/Input"
output_folder = "/home/viroicbas2023/Documents/Gmoreira/Carraças_Metagenomica_Gui/Reads_DeRNAfied/Kraken2_Experimentar"
kraken_db = "/home/viroicbas2023/Documents/Gmoreira/krakenDB/PlusPfP_GrandalhonaBro"

# Toggle exporting of unclassified reads. Potentially useful for uncharacterized sequences; Mostly useful for trash reads. True or Flase. @GmoreiraVet
EXPORT_UNCLASSIFIED = False
# Toggle memory mapping (needed when dtabase is too big for available RAM). True or False. @GmoreiraVet
USE_MEMORY_MAPPING = True

# Optional Bracken post-processing. Set RUN_BRACKEN = True to run Bracken automatically after Kraken2.
RUN_BRACKEN = True
BRACKEN_BIN = "bracken"
BRACKEN_DB = kraken_db
BRACKEN_READ_LENGTH = "1000"
BRACKEN_TAX_LEVEL = "S"
BRACKEN_THRESHOLD = "5"

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

if not os.path.isdir(input_folder):
    raise FileNotFoundError(f"Input folder not found: {input_folder}")

# Kraken2 parameters.
# Default Parameters: --threads 30 (for Watson & Crick only, check your system for available threads), --minimum-base-quality 8, --confidence 0, --memory-mapping
threads = "30"
min_base_quality = "8"
confidence = "0"
memory_mapping = "--memory-mapping"

SUPPORTED_INPUT_EXTENSIONS = [".fastq.gz", ".fq.gz", ".fastq", ".fq", ".fasta", ".fa"]

def sample_name_from_path(filename):
    for ext in sorted(SUPPORTED_INPUT_EXTENSIONS, key=len, reverse=True):
        if filename.endswith(ext):
            return filename[:-len(ext)]
    return os.path.splitext(filename)[0]

# Iterate over all FASTQ/FASTA files in the input folder
# os.listdir() lists all files in a folder.
# File name is checked for supported extensions. Only these will be affected.
input_files = sorted(
    [name for name in os.listdir(input_folder)
     if any(name.endswith(ext) for ext in SUPPORTED_INPUT_EXTENSIONS)]
)

if not input_files:
    print(f"No supported FASTQ/FASTA files found in input folder: {input_folder}")
    print("Supported extensions:", ", ".join(SUPPORTED_INPUT_EXTENSIONS))
    print("Found files:", ", ".join(sorted(os.listdir(input_folder))) or "(none)")
    raise SystemExit(1)

print(f"Found {len(input_files)} supported input file(s) in: {input_folder}\n")

for file_name in input_files:
    # Input path is constructed by joining folder path with file name.
    input_path = os.path.join(input_folder, file_name)
    # Extract the base name of the file (without extension) to use in output.
    sample_name = sample_name_from_path(file_name)
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
        "--confidence", confidence,
    ]

    # Optional export of unclassified reads. If EXPORT_UNCLASSIFIED is set to True, the script constructs an additional output path for unclassified reads and appends the necessary arguments to the Kraken2 command to export these reads to a separate FASTQ file. This allows users to retain unclassified sequences for further analysis or quality control, which can be particularly useful for identifying novel or poorly characterized organisms in metagenomic samples.
    if EXPORT_UNCLASSIFIED:
        unclassified_out = os.path.join(output_folder, f"{sample_name}_unclassified.fastq")
        cmd.extend(["--unclassified-out", unclassified_out])
    if USE_MEMORY_MAPPING:
        cmd.append("--memory-mapping")

    # Input file. The input FASTQ file is appended to the command list, which is the last argument for Kraken2. This ensures that Kraken2 processes the correct file for each iteration of the loop, allowing for efficient batch processing of multiple samples in the specified input folder.
    cmd.append(input_path)

    # Run Kraken2 command. The script prints a message indicating which file is being processed, executes the Kraken2 command using subprocess.run(), and then prints a message upon completion of each file, including the paths to the output and report files. This provides feedback to the user about the progress of the analysis and allows for easy tracking of results.
    print(f"Running Kraken2 on {file_name}...")
    subprocess.run(cmd)
    print(f"Finished {file_name}. Output: {output_path}, Report: {report_path}\n")

    if RUN_BRACKEN:
        bracken_output = os.path.join(output_folder, f"{sample_name}_bracken.tsv")
        bracken_cmd = [
            BRACKEN_BIN,
            "-d", BRACKEN_DB,
            "-i", report_path,
            "-o", bracken_output,
            "-r", BRACKEN_READ_LENGTH,
            "-l", BRACKEN_TAX_LEVEL,
        ]
        if BRACKEN_THRESHOLD:
            bracken_cmd.extend(["-t", BRACKEN_THRESHOLD])

        print(f"Running Bracken on {file_name}...")
        subprocess.run(bracken_cmd)
        print(f"Finished Bracken for {file_name}. Output: {bracken_output}\n")

print("All samples processed.")
