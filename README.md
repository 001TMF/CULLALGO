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
# NetSolP Installation (https://github.com/tvinet/NetSolP-1.0)
```bash
conda create -n CULLALGO
```
```bash
mkdir CULL
```
```bash
wget https://services.healthtech.dtu.dk/services/NetSolP-1.0/netsolp-1.0.ALL.tar.gz
```
```bash
tar -xzvf netsolp-1.0.ALL.tar.gz
```
```bash
pip install -r requirements.txt
```
```bash
#Example code. Practical implementation shown in CULLALGO1.0.py
python predict.py --FASTA_PATH ./test_fasta.fasta --OUTPUT_PATH ./test_preds.csv --MODEL_TYPE ESM12 --PREDICTION_TYPE S
```
# TemStaPro Installation (https://github.com/ievapudz/TemStaPro)
```bash
git clone https://github.com/ievapudz/TemStaPro.git
```

```bash
cd TemStaPro
```

```bash
conda env create -f environment_CPU.yml
```

```bash
conda activate temstapro_env_CPU
```

```bash
#Example usage
./temstapro -f ./tests/data/long_sequence.fasta -d ./ProtTrans/ \
    -e /home/user/CULL --mean-output ./long_sequence_predictions.tsv
```

# CULLALGO1.0
1. Have TemStaPro .tsv file output from the TemStaPro model in directory
2. CULLALGO environment (PYTHON == 3.11)
```bash
conda activate CULLALGO
```
3. Change directories in script manually
```python
# Lines 288, 289, 290
fasta_path = '/home/s_gus/progs/D.fasta'
output_path = '/home/s_gus/progs/'
thermo_path = '/home/s_gus/progs/long_sequence_predictions.tsv'

#Lines 309, 310, 311
fasta_path = '/home/s_gus/progs/D.fasta'
output_path = '/home/s_gus/progs/'    
thermo_path = "/home/s_gus/progs/long_sequence_predictions.tsv"
```
4. Run script
```bash
python CULLALGO1.0.py
```

5. After running the script, the selected_sequences.fasta file is created in the output_path directory.
```bash
nano selected_sequences.fasta
```

# Tweaks for desired results
There are multiple lines in the code that can be tweaked for individual variance. As of CULLALGO1.0, all of them are made clear via comments in the code.
