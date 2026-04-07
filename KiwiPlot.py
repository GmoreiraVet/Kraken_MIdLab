import os
import subprocess
import pysam
import matplotlib.pyplot as plt
import numpy as np
import random
import logging

# =============================================================================
# CONFIGURATION - Set all paths here.
# =============================================================================

# =============================================================================
# CONFIGURATION - Set all paths here.
# =============================================================================

OUTPUT_FOLDER = "/home/viroicbas2023/Documents/Gmoreira/EnteroVirusCarraçasUropigeaCustom/Teste"
REFERENCE_GENOME = "/home/viroicbas2023/Documents/Gmoreira/EnteroVirusCarraçasUropigeaCustom/BoleTickVirus/TIBET_BOLETICK4_REFERENCIA.fasta"
INPUT_FASTQS = [
    "/home/viroicbas2023/Documents/Gmoreira/EnteroVirusCarraçasUropigeaCustom/JXBXKX_fastq/JXBXKX_7_8T.fastq.gz",  # .gz files are supported
]

CUSTOM_COLOR = "#81B7A6"       # Set a hex color e.g. "#2A9D8F" to override, or leave as None for random

# =============================================================================
# =============================================================================

# Setup logging. This function will save logs in a file and in the console simultaneously. "logging.basicConfig" configures the logging system to write logs to a file named "pupi_log.txt" in the output folder, with a specific format that includes the timestamp, log level, and message. "logging.getLogger().addHandler(logging.StreamHandler())" adds a stream handler to the logger, which allows log messages to be printed to the console as well. "format='%(asctime)s - %(levelname)s - %(message)s'" specifies the format of the log messages, including the timestamp, log level (e.g., INFO, ERROR), and the actual log message. This setup ensures that all important information about the pipeline's execution is recorded both in a log file and displayed in real-time in the console for easy monitoring.
def setup_logging(output_folder):
    log_file = os.path.join(output_folder, "pupi_log.txt")
    logging.basicConfig(filename=log_file, level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')
    logging.getLogger().addHandler(logging.StreamHandler())

# Helper funtion. logging.info is used to log informational messages about the progress and status of the pipeline. "run_command" is a helper function that takes a command as input, logs the command being executed, and then runs it using "subprocess.run". The "shell=True" argument allows the command to be executed through the shell, and "check=True" ensures that an exception is raised if the command returns a non-zero exit status, which helps in error handling. This function centralizes the execution of shell commands and ensures that all commands are logged for better traceability and debugging.
def run_command(command):
    logging.info(f"Running command: {command}")
    subprocess.run(command, shell=True, check=True)

# This function takes a BAM file as input and calculates various coverage statistics. It uses the "pysam" library to read the BAM file and extract information about the reference genome length, total reads, mapped reads, and coverage at each position. The function calculates the average depth of coverage across the entire genome, the average depth for covered regions only, the maximum depth, and the breadth of coverage (percentage of positions covered by at least one read). It also logs all these statistics for later reference. Finally, it returns an array containing the coverage depth at each position in the reference genome.

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

    covered_positions = np.sum(coverage > 0)
    average_depth_all = np.mean(coverage)
    average_depth_covered = np.mean(coverage[coverage > 0]) if covered_positions > 0 else 0
    max_depth = np.max(coverage)
    coverage_breadth = (covered_positions / ref_length) * 100
    mapped_pct = (mapped_reads / total_reads * 100) if total_reads > 0 else 0

    logging.info(f"Reference length: {ref_length}")
    logging.info(f"Total reads: {total_reads}")
    logging.info(f"Mapped reads: {mapped_reads}")
    logging.info(f"Percentage of mapped reads: {mapped_pct:.2f}%")
    logging.info(f"Average depth (entire genome): {average_depth_all:.2f}")
    logging.info(f"Average depth (covered regions only): {average_depth_covered:.2f}")
    logging.info(f"Max depth: {max_depth}")
    logging.info(f"Coverage breadth (>=1x): {coverage_breadth:.2f}%")

    return coverage

# This function takes the raw coverage data and applies a smoothing technique using a moving average. The "np.convolve" function is used to perform the convolution of the coverage array with a window of ones, which effectively calculates the average coverage over a specified window size (default is 50). The "mode='valid'" argument ensures that only the valid part of the convolution is returned, which means that the smoothed coverage will be shorter than the original coverage by "window_size - 1" positions. This smoothed coverage can help to visualize trends in the data more clearly by reducing noise and fluctuations in the raw coverage data.

def smooth_data(coverage, window_size=50):
    return np.convolve(coverage, np.ones(window_size)/window_size, mode='valid')


# darken_color is a helper function that takes a hex color code as input and returns a darker version of that color. The function extracts the red, green, and blue components from the hex color code, multiplies each component by a specified factor (default is 0.8) to darken the color, and then constructs a new hex color code from the modified RGB values. This function is used in the plot_coverage function to create a visually appealing contrast between the line color and the fill color in the coverage plot.

def darken_color(color, factor=0.8):
    r, g, b = [int(color[i:i+2], 16) for i in (1, 3, 5)]
    r, g, b = int(r * factor), int(g * factor), int(b * factor)
    return f'#{r:02x}{g:02x}{b:02x}'

# This function generates a plot of the genome coverage depth. It takes the raw coverage and smoothed coverage as input, along with the output folder where the plot will be saved. The function randomly selects a color from a predefined list of colors for the plot, and then creates a line plot of the smoothed coverage with a filled area underneath it to visually represent the coverage depth. The plot is styled with titles, labels, grid lines, and a legend for clarity. Finally, the plot is saved as "coverage_plot.png" in the specified output folder and displayed on the screen. Logging is used to indicate where the plot has been saved for easy reference.

def plot_coverage(coverage, smoothed_coverage, output_folder):
    colors = ['#264653', '#2A9D8F', '#F4A261', '#E76F51', '#0081A7',
              '#F77F00', '#FAA307', '#FF595E', '#9A031E', '#0F4C5C']
    
    selected_color = CUSTOM_COLOR if CUSTOM_COLOR else random.choice(colors)
    line_color = darken_color(selected_color, factor=0.8)
    # ... rest of function unchanged

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

def print_ascii_kitty():
    kitty = r"""
     ∧ ,,, ∧
    ( ̳• · • ̳)
    /   づ♡ ATGCGCTGACGCAGACATAGACGACACCACACCACGGATTTAGACAGTACAGATAGGAC      


    """
    print("Kiwi has assembled your genome for you! Hard work...")
    print(kitty)

# Main function that orchestrates the entire pipeline. It starts by creating the output directory and setting up logging. It then checks if multiple FASTQ files are provided and logs the configuration details. The reference genome is indexed using minimap2, and the input FASTQ files are aligned to the reference genome. The resulting SAM file is converted to a sorted BAM file, which is then indexed for efficient access. Coverage statistics are calculated from the BAM file, and a coverage plot is generated and saved in the output folder. Finally, a message is printed to indicate that the pipeline has completed successfully, along with an ASCII art of a kitty for a fun touch!
def main():
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    setup_logging(OUTPUT_FOLDER)

    combine = len(INPUT_FASTQS) > 1
    logging.info(f"Output folder: {OUTPUT_FOLDER}")
    logging.info(f"Reference genome: {REFERENCE_GENOME}")
    logging.info(f"Input FASTQs: {INPUT_FASTQS}")
    logging.info(f"Combining FASTQs: {combine}")

    reference_mmi = os.path.join(OUTPUT_FOLDER, "reference.mmi")
    run_command(f"minimap2 -d {reference_mmi} {REFERENCE_GENOME}")

    if combine:
        combined_fastq = os.path.join(OUTPUT_FOLDER, "combined.fastq.gz")
        run_command(f"zcat {' '.join(INPUT_FASTQS)} | gzip > {combined_fastq}")  # handles both .gz and plain .fastq
        alignment_file = os.path.join(OUTPUT_FOLDER, "combined_alignment.sam")
        run_command(f"minimap2 -ax map-ont {reference_mmi} {combined_fastq} > {alignment_file}")
    else:
        alignment_file = os.path.join(OUTPUT_FOLDER, "alignment.sam")
        run_command(f"minimap2 -ax map-ont {reference_mmi} {INPUT_FASTQS[0]} > {alignment_file}")

    sorted_bam_file = os.path.join(OUTPUT_FOLDER, "alignment.sorted.bam")
    run_command(f"samtools view -Sb {alignment_file} | samtools sort -o {sorted_bam_file}")
    run_command(f"samtools index {sorted_bam_file}")

    coverage = get_coverage_and_stats(sorted_bam_file)
    smoothed_coverage = smooth_data(coverage)
    plot_coverage(coverage, smoothed_coverage, OUTPUT_FOLDER)

    logging.info(f"Pipeline completed successfully! All output files are stored in: {OUTPUT_FOLDER}")
    print_ascii_kitty()

# Entry point of the script. This ensures that the main function is called only when the script is executed directly, and not when it is imported as a module in another script. This is a common Python convention to allow for better modularity and reusability of code. (nao percebo nada disto, mas ok)
if __name__ == "__main__":
    main()