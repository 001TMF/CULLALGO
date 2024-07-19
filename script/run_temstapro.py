import subprocess
import os
from pathlib import Path
import yaml
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Process command line arguments.')
    parser.add_argument('--config', type=str, default='config.yaml', help='Path to the configuration file')
    return parser.parse_args()


# Function to read configuration from YAML file
def load_config(file_path):
    with open(file_path, 'r') as file:
        return yaml.safe_load(file)

def manage_conda_env():
    # Check the current Conda environment
    env = subprocess.run("conda info --json", shell=True, capture_output=True, text=True)
    current_env = env.stdout

    # Determine the current environment name
    if "active_env_name" in current_env:
        env_name = current_env["active_env_name"]
    else:
        # If there is no active environment, default to 'base'
        env_name = "base"

    # Manage environment based on the current active environment
    if env_name == "tempstapro_env":
        print("tempstapro_env is already activated.")
    elif env_name == "base":
        print("Activating tempstapro_env from base.")
        subprocess.run("conda activate tempstapro_env", shell=True)
    else:
        print(f"Deactivating {env_name} and activating tempstapro_env.")
        subprocess.run("conda deactivate", shell=True)
        subprocess.run("conda activate tempstapro_env", shell=True)


def run_temstapro(fasta_directory, prottrans_directory, cache_directory, output_directory):
    # Ensure the correct Conda environment is activated
    manage_conda_env()

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
            "\ -e", cache_directory,
            "--mean-output", output_path
        ]

        # Run the command and wait for it to complete
        print(f"Running command for file: {fasta_file}")
        subprocess.run(command, check=True)
        print(f"Completed processing: {fasta_file}")



def main():
    args = parse_args()
    config = load_config(args.config)
    prottrans_directory = "./ProtTrans/",
    cache_directory = "./tests/cache/",
    fasta_directory = config["TSP_paths"]["fasta_directory"]
    output_directory = config["TSP_paths"]["output_directory"]
    run_temstapro(fasta_directory, prottrans_directory, cache_directory, output_directory)

if __name__ == "__main__":
    main()
