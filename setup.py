import subprocess
import os

def run_command(command, cwd=None, capture_output=False):
    """Utility function to run shell commands with optional output capture and error handling."""
    try:
        if capture_output:
            result = subprocess.run(command, check=True, shell=True, cwd=cwd, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, text=True)
            return True, result.stdout
        else:
            subprocess.run(command, check=True, shell=True, cwd=cwd)
            return True, None
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
        return False, e.output if capture_output else None


def is_conda_environment_present(env_name):
    """Check if a conda environment is already present by parsing `conda env list`."""
    success, output = run_command("conda env list", capture_output=True)
    if success and output:
        return env_name in output
    else:
        print("Could not list conda environments.")
        return False


def clone_repo(url, directory):
    """Clone a git repository if the directory does not exist."""
    if not os.path.exists(directory):
        print(f"Cloning {url} into {directory}...")
        return run_command(f"git clone {url} {directory}")
    else:
        print(f"Directory {directory} already exists. Skipping clone.")
        return False


def setup_conda_environment(env_name, requirements_file):
    """Setup a conda environment if it's not already present."""
    if not is_conda_environment_present(env_name):
        print(f"Creating Conda environment {env_name}...")
        if run_command(f"conda create -n {env_name} python=3.11 -y"):
            print(f"Installing requirements from {requirements_file}...")
            run_command(f"pip install -r {requirements_file}")
    else:
        print(f"Conda environment {env_name} already exists. Skipping creation.")


def main():
    print("Setting up the CULLALGO environment...")

    setup_conda_environment("CULLALGO", "cull_requirements.txt")

    if input("Would you like to install NETSOLP? (yes/no): ").strip().lower() == 'yes':
        if not os.path.exists("netsolp-1.0.ALL.tar.gz"):
            run_command("wget https://services.healthtech.dtu.dk/services/NetSolP-1.0/netsolp-1.0.ALL.tar.gz")
        if run_command("tar -xzvf netsolp-1.0.ALL.tar.gz"):
            if os.path.exists("requirements.txt"):
                run_command("pip install -r requirements.txt")
            else:
                print("No requirements file found. Manual installation will be required")
    else:
        print("Skipping installation of NETSOLP.")

    if input("Would you like to install TemStaPro? (yes/no): ").strip().lower() == 'yes':
        temsta_directory = "./temstapro"
        if clone_repo("https://github.com/ievapudz/TemStaPro.git", temsta_directory):
            os.chdir(temsta_directory)
            setup_type = input("Do you want to set up for GPU or CPU? Enter 'GPU' or 'CPU': ").strip().upper()
            env_file = "environment_GPU.yml" if setup_type == "GPU" else "environment_CPU.yml"
            env_name = "temstapro_env_GPU" if setup_type == "GPU" else "temstapro_env_CPU"
            if not is_conda_environment_present(env_name):
                run_command(f"conda env create -f {env_file}")
            else:
                print(f"Environment {env_name} already exists. Skipping setup.")
        else:
            print("Error or user cancelled cloning TemStaPro. Cannot proceed with its setup.")
    else:
        print("Skipping installation of TemStaPro.")

    print("Setup completed. Please activate the conda environments as needed to use the tools.")


if __name__ == "__main__":
    main()
