import subprocess
import os
from pathlib import Path

# -----------------------------
# CONFIGURATION (edit here)
# -----------------------------

FASTA_DIR = "/home/viroicbas2023/Documents/Gmoreira/Vibrios/VibrioCHolerae/New/"
OUTPUT_DIR = "/home/viroicbas2023/Documents/Gmoreira/Vibrios/VibrioCHolerae/New/ABRICATE"

# -----------------------------
# SCRIPT START
# -----------------------------

# Make sure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Get list of installed Abricate databases
result = subprocess.run(["abricate", "--list"], capture_output=True, text=True)
lines = result.stdout.strip().split("\n")

# Skip the header line
databases = [line.split()[0] for line in lines[1:]]

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