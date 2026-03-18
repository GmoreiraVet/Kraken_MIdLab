from Bio import SeqIO
import os
import re

# =========================
# Variables - Update before Running
# =========================
INPUT_FASTA = "/home/viroicbas2023/Documents/Gmoreira/Vibrios/Ecoli/NEW/Virsorter2/final-viral-combined.fa"
OUTPUT_DIR = "/home/viroicbas2023/Documents/Gmoreira/Vibrios/Ecoli/NEW/Virsorter2/FinalViral-separated"
USE_FULL_HEADER = True  # True = use full header, False = only first word

 # Function to sanitize filenames. Replaces any nonalphanumeric character with underscores. Avoids compatibility issues.
def sanitize_filename(name):
    """Make safe filenames from FASTA headers"""
    return re.sub(r'[^A-Za-z0-9_.-]', '_', name)

# Splits Multifasta into individual fasta files. Each file is named after the header of the sequence, sanitized to ensure it is a valid filename. If USE_FULL_HEADER is set to True, the entire header is used for naming; otherwise, only the first word of the header is used. The function iterates through each record in the input FASTA file, creates an output file for each sequence, and writes the sequence to its respective file in the specified output directory.
def split_fasta(input_fasta, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for i, record in enumerate(SeqIO.parse(input_fasta, "fasta"), start=1):
        
        if USE_FULL_HEADER:
            name = record.description
        else:
            name = record.id

        name = sanitize_filename(name)

        # fallback if header is weird/empty
        if not name:
            name = f"contig_{i}"

        output_file = os.path.join(output_dir, f"{name}.fasta")

        with open(output_file, "w") as out_f:
            SeqIO.write(record, out_f, "fasta")

    print(f"✅ Done. Files saved in: {output_dir}")


# Run
split_fasta(INPUT_FASTA, OUTPUT_DIR)