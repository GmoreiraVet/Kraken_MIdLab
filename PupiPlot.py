import os
import subprocess
import pysam
import matplotlib.pyplot as plt
import numpy as np
import random
import logging

# Set up logging
def setup_logging(output_folder):
    log_file = os.path.join(output_folder, "pupi_log.txt")
    logging.basicConfig(filename=log_file, level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    logging.getLogger().addHandler(logging.StreamHandler())

def run_command(command):
    logging.info(f"Running command: {command}")
    subprocess.run(command, shell=True, check=True)

def get_coverage_and_stats(bam_file):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    ref_length = samfile.lengths[0]
    coverage = np.zeros(ref_length, dtype=int)

    total_reads = 0
    mapped_reads = 0

    for read in samfile.fetch(until_eof=True):
        total_reads += 1
        if not read.is_unmapped:
            mapped_reads += 1

    for pileupcolumn in samfile.pileup():
        pos = pileupcolumn.reference_pos
        coverage[pos] = pileupcolumn.n

    samfile.close()

    # Statistics
    covered_positions = np.sum(coverage > 0)
    average_depth_all = np.mean(coverage)
    average_depth_covered = np.mean(coverage[coverage > 0]) if covered_positions > 0 else 0
    max_depth = np.max(coverage)
    coverage_breadth = (covered_positions / ref_length) * 100
    mapped_pct = (mapped_reads / total_reads * 100) if total_reads > 0 else 0

    # Logging
    logging.info(f"Reference length: {ref_length}")
    logging.info(f"Total reads: {total_reads}")
    logging.info(f"Mapped reads: {mapped_reads}")
    logging.info(f"Percentage of mapped reads: {mapped_pct:.2f}%")
    logging.info(f"Average depth (entire genome): {average_depth_all:.2f}")
    logging.info(f"Average depth (covered regions only): {average_depth_covered:.2f}")
    logging.info(f"Max depth: {max_depth}")
    logging.info(f"Coverage breadth (>=1x): {coverage_breadth:.2f}%")

    return coverage

def smooth_data(coverage, window_size=50):
    return np.convolve(coverage, np.ones(window_size)/window_size, mode='valid')

def darken_color(color, factor=0.8):
    r, g, b = [int(color[i:i+2], 16) for i in (1, 3, 5)]
    r, g, b = int(r * factor), int(g * factor), int(b * factor)
    return f'#{r:02x}{g:02x}{b:02x}'

def plot_coverage(coverage, smoothed_coverage, output_folder):
    colors = ['#264653', '#2A9D8F', '#F4A261', '#E76F51', '#0081A7',
              '#F77F00', '#FAA307', '#FF595E', '#9A031E', '#0F4C5C']
    selected_color = random.choice(colors)
    line_color = darken_color(selected_color, factor=0.8)

    plt.figure(figsize=(10, 6))
    plt.plot(smoothed_coverage, color=line_color, label='Smoothed Coverage')
    plt.fill_between(range(len(smoothed_coverage)), smoothed_coverage,
                     color=selected_color, alpha=0.6)
    plt.title("Genome Coverage Depth", color='#333333')
    plt.xlabel("Position on Genome", color='#333333')
    plt.ylabel("Coverage Depth", color='#333333')
    plt.grid(True, color='#333333', alpha=0.3)
    plt.legend()

    output_path = os.path.join(output_folder, "coverage_plot.png")
    plt.savefig(output_path)
    plt.show()
    logging.info(f"Coverage plot saved to: {output_path}")

def print_ascii_bunny():
    bunny = r"""
      /\ /|
      \ V/
      | "")
      /  |
     /  \
    *(__\_\)
    """
    print("Pupi has assembled your genome for you!")
    print(bunny)

def main():
    output_folder = input("Enter the name of the output folder: ").strip()
    os.makedirs(output_folder, exist_ok=True)
    setup_logging(output_folder)

    reference_genome = input("Enter the path to the reference genome (reference.fasta): ").strip()
    input_fastqs = input("Enter the path(s) to the input FASTQ file(s) (comma separated if more than one): ").strip().split(',')
    combine_fastqs = input("Do you want to combine the FASTQ files? (yes/no): ").strip().lower()

    reference_mmi = os.path.join(output_folder, "reference.mmi")
    run_command(f"minimap2 -d {reference_mmi} {reference_genome}")

    combined_fastq = os.path.join(output_folder, "combined.fastq")
    if combine_fastqs == "yes":
        run_command(f"cat {' '.join([fastq.strip() for fastq in input_fastqs])} > {combined_fastq}")

    if combine_fastqs == "yes":
        alignment_file = os.path.join(output_folder, "combined_alignment.sam")
        run_command(f"minimap2 -ax map-ont {reference_mmi} {combined_fastq} > {alignment_file}")
    else:
        alignment_file = os.path.join(output_folder, "alignment.sam")
        run_command(f"minimap2 -ax map-ont {reference_mmi} {input_fastqs[0].strip()} > {alignment_file}")

    sorted_bam_file = os.path.join(output_folder, "alignment.sorted.bam")
    run_command(f"samtools view -Sb {alignment_file} | samtools sort -o {sorted_bam_file}")
    run_command(f"samtools index {sorted_bam_file}")

    coverage = get_coverage_and_stats(sorted_bam_file)
    smoothed_coverage = smooth_data(coverage)
    plot_coverage(coverage, smoothed_coverage, output_folder)

    logging.info(f"Pipeline completed successfully! All output files are stored in: {output_folder}")
    print_ascii_bunny()

if __name__ == "__main__":
    main()
