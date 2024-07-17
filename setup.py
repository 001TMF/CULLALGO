import subprocess
import os

def run_command(command):
    """Utility function to run shell commands."""
    try:
        subprocess.run(command, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")

def confirm_and_clone(url, directory):
    """Confirm before cloning a git repository."""
    response = input(f"Do you want to download and set up from {url}? (yes/no): ").strip().lower()
    if response == 'yes':
        run_command(f"git clone {url} {directory}")


def main():
    print("Setting up the CULLALGO environment...")

    # Setting up CULLALGO environment
    print("Creating Conda environment for CULLALGO...")
    run_command("conda create -n CULLALGO python=3.11 -y")
    run_command("conda activate CULLALGO")
    print("Installing CULLALGO requirements...")
    run_command("pip install -r cull_requirements.txt")

    # Install NETSOLP within the cullalgo directory
    nsp_install = print("would you like to install netsolp? (yes/no): ")
    if nsp_install == 'yes':
        print("Downloading and setting up NETSOLP...")
        run_command("mkdir netsolp")
        os.chdir("netsolp")
        run_command("wget https://services.healthtech.dtu.dk/services/NetSolP-1.0/netsolp-1.0.ALL.tar.gz")
        run_command("tar -xzvf netsolp-1.0.ALL.tar.gz")
        run_command("pip install -r requirements.txt")
        os.chdir("..")
    else:
        print("Skipping installation of netsolp")

    # Setting up TemStaPro
    print("Setting up TemStaPro...")
    temsta_directory = "./temstapro"  # Specify TemStaPro directory
    confirm_and_clone("https://github.com/ievapudz/TemStaPro.git", temsta_directory)
    os.chdir(temsta_directory)

    # Ask user for setup choice (GPU or CPU)
    setup_type = input("Do you want to set up for GPU or CPU? Enter 'GPU' or 'CPU': ").strip().upper()
    if setup_type == "GPU":
        run_command("conda deactivate")
        run_command("conda env create -f environment_GPU.yml")
        run_command("conda activate temstapro_env_GPU")
    else:
        run_command("conda deactivate")
        run_command("conda env create -f environment_CPU.yml")
        run_command("conda activate temstapro_env_CPU")

    print("Setup completed. Please activate the conda environments as needed to use the tools.")

if __name__ == "__main__":
    main()
