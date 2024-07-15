# CULLALGO1.0
'CULLALGO is a script that uses multiple parameters to 'cull' a large dataset of .fasta protein sequences to a manageable amount with desired traits in a .fasta output. All measurements are based on up to date existing literature.'

The determinant parameters in CULLALGO1.0 are as follows:
1. Molecular weight
2. Average Surface Accesibility
3. Isoelectric Point
4. Cost
5. Sollubility
6. Thermostability
7. Alpha helical propensity
8. Shannon's entropy (Redundancy measure)
9. DNA complexity

'CULLALGO uses NETSOLP and TemStaPro for determining Sollubility and Thermostability thresholds, respectively. Both are required to install and set up in order to run the script.'

'NETSOLP https://github.com/tvinet/NetSolP-1.0'

'TemStaPro https://github.com/ievapudz/TemStaPro'

# Set-up
# NetSolP Installation from (https://github.com/tvinet/NetSolP-1.0)
```bash
conda create -n CULLALGO
```
```bash
mkdir CULL
```
```bash
wget https://services.healthtech.dtu.dk/services/NetSolP-1.0/netsolp-1.0.ALL.tar.gz
```
1. Set up a seperate conda environment with all required python imports from the requirements.txt file. Make sure your python version is AT LEAST 3.11. If requirements.txt poses problems, run script manually and /pip install/ whatever is needed. (WARNING: A seperate environment, where python==3.7 needs to exist to run TemStaPro. This is mentioned in its installation guidelines.)
2. Install and correctly set-up NETSOLP and TemStaPro. It is advised to first test both programs independently to make sure they work. As of CULLALGO1.0, NETSOLP Sollubility measurements CAN be run natively in the script after set-up. However, TemStaPro MUST BE FIRST RUN SEPERATELY, and the output file correctly placed in the desired directory. TemStaPro cannot be run natively in the script.
3. All directories must be changed manually to match desired output. For error handling, all files used for CULLALGO1.0 should be in the same directory.
![image](https://github.com/user-attachments/assets/be323027-2965-4be6-a8aa-5de9d1994005)
![image](https://github.com/user-attachments/assets/2445f9c4-1cf6-4af9-bda1-12b36c36ac22)
These are the lines that need to be changed before running the script, in the same format.
4. After running the script, the selected_sequences.fasta file is created in the output_path directory.

# Tweaks for desired results
There are multiple lines in the code that can be tweaked for individual variance. As of CULLALGO1.0, all of them are made clear via comments in the code.
