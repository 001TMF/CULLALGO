import os
import subprocess

def run_command(command, cwd=None, capture_output=False):
    """Utility function to run shell commands with optional output capture and error handling."""
    if capture_output:
        try:
            if capture_output:
                result = subprocess.run(command, check=True, shell=True, cwd=cwd, stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE, text=True)
                return True, result.stdout
            else:
                subprocess.run(command, check=True, shell=True, cwd=cwd)
                return True, None
        except Exception as e:
            print(f"An error occurred: {e}")
            return False, str(e)
    else:
        # If show_progress is True, do not capture stdout so it prints directly to the console.
        try:
            process = subprocess.Popen(command, shell=True, cwd=cwd)
            process.wait()  # Wait for the process to complete
            return process.returncode == 0, None
        except Exception as e:
            print(f"An error occurred: {e}")
            return False, str(e)

def clone_repo(url, directory):
    """Clone a git repository if the directory does not exist."""
    if not os.path.exists(directory):
        print(f"Cloning {url} into {directory}...")
        success, output = run_command(f"git clone {url} {directory}")
        if not success:
            print(output)
        return success
    else:
        print(f"Directory {directory} already exists. Skipping clone.")
        return False

def setup_conda_environment(env_name, requirements_file):
    """Setup a conda environment and install requirements."""
    print(f"Checking for Conda environment {env_name}...")
    success, output = run_command("conda env list", capture_output=True)
    print(output)
    if env_name not in output:
        print(f"Creating Conda environment {env_name}...")
        run_command(f"conda create -n {env_name} python=3.11 -y")
    print(f"{env_name} already exists. refreshing requirements")
    print(f"Installing requirements from {requirements_file} in {env_name}...")
    run_command(f"conda run -n {env_name} pip install -r {requirements_file}")

def main():
    print("Setting up the CULLALGO environment...")
    setup_conda_environment("CULLALGO", "cull_requirements.txt")

    if input("Would you like to install NETSOLP? (yes/no): ").strip().lower() == 'yes':
        if not os.path.exists("netsolp"):
            os.makedirs("netsolp")
        if not os.path.exists("netsolp/netsolp-1.0.ALL.tar.gz"):
            # Run wget without output capture to display progress
            success, output = run_command("wget https://services.healthtech.dtu.dk/services/NetSolP-1.0/netsolp-1.0.ALL.tar.gz")
            if not success:
                print(output or "Failed to download NETSOLP.")
        if os.path.exists("netsolp/netsolp-1.0.ALL.tar.gz"):
            os.chdir("netsolp")
            success, output = run_command("tar -xzvf netsolp-1.0.ALL.tar.gz")
            os.chdir('..')
            if not success:
                print(output or "Failed to extract NETSOLP.")
            if os.path.exists("netsolp/requirements.txt"):
                run_command(f"conda run -n CULLALGO pip install -r requirements.txt")
                print("Installing requirements from requirements.txt")
            else:
                print("No requirements file found. Manual installation may be required.")
    else:
        print("Skipping installation of NETSOLP.")

    if input("Would you like to install TemStaPro? (yes/no): ").strip().lower() == 'yes':
        temsta_directory = "./temstapro"
        if clone_repo("https://github.com/ievapudz/TemStaPro.git", temsta_directory):
            setup_type = input("Do you want to set up for GPU or CPU? Enter 'GPU' or 'CPU': ").strip().upper()
            env_file = "environment_GPU.yml" if setup_type == "GPU" else "environment_CPU.yml"
            env_name = "temstapro_env_GPU" if setup_type == "GPU" else "temstapro_env_CPU"
            run_command(f"conda env create -f {temsta_directory}/{env_file}")
        else:
            print("Error or user cancelled cloning TemStaPro. Cannot proceed with its setup.")
    else:
        print("Skipping installation of TemStaPro.")

    print("Setup completed. Please use the conda environments as needed to use the tools.")

if __name__ == "__main__":
    main()
