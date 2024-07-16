# CULLALGO 1.1

## Description
CULLALGO is a script designed to "cull" large datasets of .fasta protein sequences into a manageable amount with desired traits, outputting results in a .fasta file. All measurements leverage up-to-date literature.

### Determinant Parameters
The script uses the following parameters:
1. Molecular Weight
2. Average Surface Accessibility
3. Isoelectric Point
4. Cost
5. Solubility
6. Thermostability
7. Alpha Helical Propensity
8. Shannon's Entropy (Redundancy Measure)
9. DNA Complexity

CULLALGO utilizes both NETSOLP and TemStaPro for determining solubility and thermostability thresholds, respectively. Both tools must be installed and configured to run the script.

### Dependencies
- NETSOLP: [GitHub - NetSolP](https://github.com/tvinet/NetSolP-1.0)
- TemStaPro: [GitHub - TemStaPro](https://github.com/ievapudz/TemStaPro)

## Automatic installation
### Important Notes
#### ⚠️ Still in development ⚠️
- **Permissions**: Running such a script may require administrative privileges, especially for parts that involve installing system-wide packages or modifying system paths.
- **Compatibility**: This script assumes a Unix-like operating system because of its dependency on bash commands. For Windows, adjustments might be necessary, particularly in how environments are activated and paths are handled.
### Run
```bash
python3 setup.py
```


## Manaual Installation
### Environment Setup
```bash
conda create -n CULLALGO python=3.11
```
```bash
conda activate CULLALGO
```

### CULLALGO Script preperation
```bash
git clone https://github.com/gusgrazelis/CULLALGO.git
```
```bash
cd cullalgo
```

### NETSOLP Installation
This should be installed within the cullalgo directory
```bash
wget https://services.healthtech.dtu.dk/services/NetSolP-1.0/netsolp-1.0.ALL.tar.gz
```
```bash
tar -xzvf netsolp-1.0.ALL.tar.gz
```
```bash
pip install -r requirements.txt
```

### Running CULLALGO
```bash
python3 CULLALGO1.1.py --config config.yaml
```

### TemStaPro Installation
This should be its own directory and can be where you like.
```bash
git clone https://github.com/ievapudz/TemStaPro.git
```
```bash
cd TemStaPro
```
#### Environment requirements

Before starting up Anaconda or Miniconda should be installed
in the system. Follow instructions given in 
[Conda's documentation.](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)

Setting up the environment can be done in one of the following ways.

#### From YML file

In this repository two YML files can be found: one YML file
has the prerequisites for the environment that exploits only 
CPU ([`environment_CPU.yml`](./environment_CPU.yml)), another one to exploit both CPU 
GPU ([`environment_GPU.yml`](./environment_GPU.yml)).

This approach was tested with Conda 4.10.3 and 4.12.0 versions.

Run the following command to create the environment from a 
YML file:
```
conda env create -f environment_CPU.yml
```

Activate the environment:
```
conda activate temstapro_env_CPU
```

#### From scratch
#### GPU Setup
To set up the environment to exploit **GPU** for the program, run the following commands:
```bash
conda create -n temstapro_env python=3.7
```
```bash
conda activate temstapro_env
```
```bash
conda install pytorch torchvision torchaudio pytorch-cuda=11.7 -c pytorch -c nvidia
```
```bash
conda install -c conda-forge transformers
```
```bash
conda install -c conda-forge sentencepiece
```
```bash
conda install -c conda-forge matplotlib
```

To test if PyTorch package is installed to exploit CUDA,
call `python3` command interpreter and run the 
following lines:
```bash
python3
```
```python
import torch
torch.cuda.is_available()
```

If the output is 'True', then the installing procedure was successful,
otherwise try to set the path to the installed packages:
```
export PATH=/usr/local/cuda-11.7/bin${PATH:+:${PATH}}
export LD_LIBRARY_PATH=/usr/local/cuda-11.7/lib64\${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
```

If CUDA for PyTorch is still not available, check out the [forum.](https://github.com/pytorch/pytorch/issues/30664)

#### CPU setup
For the systems without GPU, run the following commands for ***CPU*** setup:
```
conda create -n temstapro_env python=3.7
conda activate temstapro_env
conda install -c conda-forge transformers
conda install pytorch -c pytorch
conda install -c conda-forge sentencepiece
conda install -c conda-forge matplotlib
```