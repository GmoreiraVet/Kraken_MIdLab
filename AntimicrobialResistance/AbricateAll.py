import subprocess
import os
from pathlib import Path

# -----------------------------
# Configuration - Update paths as needed
# -----------------------------

FASTA_DIR = "/home/viroicbas2023/Documents/Gmoreira/Vibrios/VibrioCHolerae/New/"
OUTPUT_DIR = "/home/viroicbas2023/Documents/Gmoreira/Vibrios/VibrioCHolerae/New/ABRICATE"

# -----------------------------
# SCRIPT START
# -----------------------------

# Make sure output directory exists. If it doesn't, create it. This is where all the Abricate results will be stored, organized by database.
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Gets a list of all availkable abricate databases using "abricate --list". Output is parsed to extract the database names, Which are then stored in a list for later use
result = subprocess.run(["abricate", "--list"], capture_output=True, text=True)
# lines=result.stdout.strip().split("\n"). Splits the output of the command into individual lines, creating a list where each element is a line from the output. The first line is typically a header, so it is skipped in the next step. Each subsequent line is split into columns (using whitespace as a delimiter), and the first column (which contains the database name) is extracted and stored in the "databases" list. This allows the script to know which databases are available for running Abricate against the FASTA files.
lines = result.stdout.strip().split("\n")

# Skip the header line
databases = [line.split()[0] for line in lines[1:]]

# Prints found databases.
print(f"Found databases: {databases}")

# Loop through databases
for db in databases:
    print(f"\nRunning Abricate on database: {db}")
    db_output_dir = Path(OUTPUT_DIR) / db
    db_output_dir.mkdir(exist_ok=True)
    
    # Loop through FASTA files
    fasta_files = list(Path(FASTA_DIR).glob("*.fasta")) + list(Path(FASTA_DIR).glob("*.fa"))
    if not fasta_files:
        print(f"No FASTA files found in {FASTA_DIR}")
        continue
    
    for fasta_file in fasta_files:
        output_file = db_output_dir / f"{fasta_file.stem}.txt"
        print(f"  Processing {fasta_file.name} -> {output_file.name}")
        subprocess.run(["abricate", "--db", db, str(fasta_file)], stdout=open(output_file, "w"))

print("\nAll Abricate runs completed!")