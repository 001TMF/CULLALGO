import subprocess
import os
from pathlib import Path

def run_temstapro(fasta_directory, prottrans_directory, cache_directory, output_directory):
    # List all FASTA files in the given directory
    fasta_files = [file for file in os.listdir(fasta_directory) if file.endswith('.fasta')]

    for fasta_file in fasta_files:
        # Construct the full path to the FASTA file
        fasta_path = os.path.join(fasta_directory, fasta_file)

        # Generate output file path
        output_file_name = f"{Path(fasta_file).stem}-TemStaPro.tsv"
        output_path = os.path.join(output_directory, output_file_name)

        # Build the command
        command = [
            "./temstapro",
            "-f", fasta_path,
            "-d", prottrans_directory,
            "-e", cache_directory,
            "--mean-output", output_path
        ]

        # Run the command and wait for it to complete
        print(f"Running command for file: {fasta_file}")
        subprocess.run(command, check=True)
        print(f"Completed processing: {fasta_file}")

# Example usage
run_temstapro(
    fasta_directory="/home/s_tristan/externals/CULLALGO/v2_test/fasta",
    prottrans_directory="./ProtTrans/",
    cache_directory="./tests/cache/",
    output_directory="/home/s_tristan/externals/CULLALGO/v2_test/TSP"
)
