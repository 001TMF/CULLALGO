import subprocess
import json
import os

# Paths to the scripts and config file
script_dir = "./scripts"
config_file = "config.yaml"
run_temstapro_script = os.path.join(script_dir, "run_temstapro.py")
cullalgo_script = os.path.join(script_dir, "CULLALGOv2.py")


def run_command(command):
    result = subprocess.run(command, shell=True)
    return result.returncode


def get_current_conda_env():
    conda_info = subprocess.check_output("conda info --json", shell=True)
    conda_info_json = json.loads(conda_info)
    return os.path.basename(conda_info_json["active_prefix"])


def activate_conda_env(env_name):
    subprocess.run(f"conda activate {env_name}", shell=True, check=True)


def deactivate_conda_env():
    subprocess.run("conda deactivate", shell=True, check=True)


# Run the run_temstapro script with the config file
if run_command(f"{run_temstapro_script} --config {config_file}") == 0:
    print("run_temstapro completed successfully.")

    current_env = get_current_conda_env()

    if current_env not in ["base", "CULLALGO"]:
        print(f"Deactivating current environment: {current_env}")
        deactivate_conda_env()

    print("Activating CULLALGO environment.")
    activate_conda_env("CULLALGO")

    # Run the CULLALGOv2.py script with the config file
    run_command(f"python {cullalgo_script} --config {config_file}")

    # Optionally, you can deactivate the environment if needed
    # print("Deactivating CULLALGO environment.")
    # deactivate_conda_env()
else:
    print("run_temstapro failed. Exiting script without running CULLALGOv2.py.")
