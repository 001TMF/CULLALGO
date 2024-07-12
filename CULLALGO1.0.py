import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import os
import subprocess
from collections import Counter
import math

# Define weights, ASAs, costs, Alpha-helical propensities, Codon dict.
weights = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
           'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
           'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
           'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06}

ASAs = {'A': 0.405, 'C': 0.268, 'D': 0.615, 'E': 0.586, 'F': 0.290,
        'G': 0.588, 'H': 0.425, 'I': 0.273, 'K': 0.607, 'L': 0.321,
        'M': 0.364, 'N': 0.568, 'P': 0.502, 'Q': 0.573, 'R': 0.539,
        'S': 0.568, 'T': 0.480, 'V': 0.306, 'W': 0.279, 'Y': 0.319}

costs = {'A': 11.7, 'C': 24.7, 'D': 12.7, 'E': 15.3, 'F': 52.0,
         'G': 11.7, 'H': 38.3, 'I': 32.3, 'K': 30.3, 'L': 27.3,
         'M': 34.3, 'N': 14.7, 'P': 20.3, 'Q': 16.3, 'R': 27.3,
         'S': 11.7, 'T': 18.7, 'V': 23.3, 'W': 74.3, 'Y': 50.0}

alpha_helix_propensity = {
    'A': 1.45, 'R': 0.79, 'N': 0.73, 'D': 0.98,
    'C': 0.77, 'Q': 1.17, 'E': 1.53, 'G': 0.53,
    'H': 1.00, 'I': 1.00, 'L': 1.34, 'K': 1.07,
    'M': 1.20, 'F': 1.12, 'P': 0.59, 'S': 0.79,
    'T': 0.82, 'W': 1.14, 'Y': 0.61, 'V': 1.14
}

codon_dict = {
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'N': ['AAU', 'AAC'],
    'D': ['GAU', 'GAC'],
    'C': ['UGU', 'UGC'],
    'Q': ['CAA', 'CAG'],
    'E': ['GAA', 'GAG'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    'H': ['CAU', 'CAC'],
    'I': ['AUU', 'AUC', 'AUA'],
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'K': ['AAA', 'AAG'],
    'M': ['AUG'],
    'F': ['UUU', 'UUC'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'W': ['UGG'],
    'Y': ['UAU', 'UAC'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG']
}

# DNA complexity
def calculate_dna_complexity(sequence):
    total_codons = sum(len(codon_dict[aa]) for aa in sequence if aa in codon_dict)
    return total_codons

# Molecular weight
def molweight(sequence):
    return sum(weights[aa] for aa in sequence if aa in weights)

# Average Surface Accessibility
def ASA(sequence):
    sums = sum(ASAs[aa] for aa in sequence if aa in ASAs)
    TASA = sums / len(sequence) if len(sequence) > 0 else 0.0
    return TASA

# Isoelectric Point
def ISO_P(sequence):
    protein = IP(sequence)
    return protein.pi()

# Cost Parameter
def costing(sequence):
    return sum(costs[aa] for aa in sequence if aa in costs)

# Calculates average propensity scores of amino acids in a specified 'window size' of the entire sequence.
# Adjust window size if needed (default = 10)
# Adjust threshold for narrower picking (default = 1.0)
def calculate_propensity(sequence, window_size=10, threshold=1.0):
    propensities = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        avg_propensity = sum(alpha_helix_propensity.get(aa, 0) for aa in window) / window_size
        propensities.append(avg_propensity)
    return propensities

# Iterates propensities to determine helical regions
def detect_alpha_helices(sequence, window_size=10, threshold=1.0):
    propensities = calculate_propensity(sequence, window_size, threshold)
    alpha_helices = []
    start = None

    for i, propensity in enumerate(propensities):
        if propensity >= threshold:
            if start is None:
                start = i
        else:
            if start is not None:
                end = i + window_size - 1
                alpha_helices.append((start, end))
                start = None

    if start is not None:
        end = len(sequence)
        alpha_helices.append((start, end))

    return alpha_helices

def read_fasta(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append((record.id, str(record.seq)))
    return sequences

def read_solubility(file_path):
    df = pd.read_csv(file_path, header=None)
    df.columns = ['ID', 'Sequence', 'Solubility', 'Extra1', 'Extra2', 'Extra3', 'Extra4', 'Extra5']
    df['Solubility'] = pd.to_numeric(df['Solubility'], errors='coerce')
    return df[['ID', 'Solubility']]

def read_thermostability(output_path, thermo_path):
    thermo_path = os.path.join(output_path, "long_sequence_predictions.tsv")
    df = pd.read_csv(thermo_path, sep='\t')
    df.columns = [f'Col{i+1}' for i in range(len(df.columns))]
    thermostability_column = df.columns[16] 
    df['Thermostability'] = df[thermostability_column].apply(extract_last_number)
    return df[['Col1', 'Thermostability']].rename(columns={'Col1': 'ID'})
    
def extract_last_number(value):
    import numpy as np
    if isinstance(value, float) and np.isnan(value):
        return None  
    elif isinstance(value, (float, int)):
        return int(value)

    import re
    numbers = re.findall(r'\d+', str(value))
    if numbers:
        # Return the last sequence of digits found
        return int(numbers[-1])
    return None

    
# Shannon entropy (redundancy measure) (Bigger number, more complex -> Good.)
def calculate_entropy(sequence):
    length = len(sequence)
    counts = Counter(sequence)
    entropy = -sum((count/length) * math.log2(count/length) for count in counts.values())
    return entropy

# Calculates all properties in format
def calculate_properties(fasta_path, solubility_path, thermo_path):
    sequences = read_fasta(fasta_path)
    solubility_data = read_solubility(solubility_path)
    thermostability_data = read_thermostability(output_path, thermo_path)

    data = []
    for seq_id, sequence in sequences:
        weight = molweight(sequence)
        asa = ASA(sequence)
        isoelectric_point = ISO_P(sequence)
        cost = costing(sequence)
        alpha_helices = detect_alpha_helices(sequence)
        num_alpha_helices = len(alpha_helices)
        entropy = calculate_entropy(sequence)
        dna_complexity = calculate_dna_complexity(sequence)
        solubility_row = solubility_data[solubility_data['ID'] == seq_id]
        thermostability_row = thermostability_data[thermostability_data['ID'] == seq_id]

        if not solubility_row.empty:
            solubility = solubility_row['Solubility'].values[0]
            thermostability = thermostability_row['Thermostability'].values[0] if not thermostability_row.empty else None

        data.append((seq_id, weight, asa, isoelectric_point, cost, solubility, thermostability, num_alpha_helices, entropy, dna_complexity))
    
    df = pd.DataFrame(data, columns=['ID', 'Molecular_Weight', 'ASA', 'Isoelectric_Point', 'Cost', 'Solubility', 'Thermostability', 'Num_Alpha_Helices', 'Entropy', 'DNA_Complexity'])
    return df

# Threshold percentages for manual tweaking of starting thresholds
def determine_thresholds_percentile(df, metric, percentile=100, reverse=False):
    metric_values = df[metric]
    if reverse:
        min_threshold = np.percentile(metric_values, 100 - percentile)
        max_threshold = max(metric_values)
    else:
        min_threshold = min(metric_values)
        max_threshold = np.percentile(metric_values, 100 - percentile)
    return min_threshold, max_threshold

# User input prompt for percentile determination on 'general' cues. I.e., how much something matters.
def ask_user_for_weights():
    parameters = ["Molecular_Weight", "ASA", "Isoelectric_Point", "Cost", "Solubility", "Thermostability", "Num_Alpha_Helices", "Entropy", "DNA_Complexity"]
    weightings = {}
    for parameter in parameters:
        while True:
            try:
                weight = int(input(f"How much does {parameter} matter: 1 - Doesn't matter, 2 - Matters somewhat, 3 - Paramount: "))
                if weight in [1, 2, 3]:
                    if weight == 1:
                        weightings[parameter] = 0
                    elif weight == 2:
                        weightings[parameter] = 50
                    else:
                        weightings[parameter] = 80
                    break
                else:
                    print("Please enter 1, 2, or 3.")
            except ValueError:
                print("Invalid input. Please enter a number (1, 2, or 3).")
    return weightings
    
def write_fasta(sequences, output_path):
    fasta_records = []
    for seq_id, sequence in sequences:
        record = SeqRecord(Seq(sequence), id=seq_id, description="")
        fasta_records.append(record)
    SeqIO.write(fasta_records, output_path, "fasta")

#Select sequences from user thresholds
def select_sequences(df, thresholds, num_sequences):
    good_sequences = df[
        (df["Molecular_Weight"].between(thresholds["Molecular_Weight"][0], thresholds["Molecular_Weight"][1])) &
        (df["ASA"].between(thresholds["ASA"][0], thresholds["ASA"][1])) &
        (df["Isoelectric_Point"].between(thresholds["Isoelectric_Point"][0], thresholds["Isoelectric_Point"][1])) &
        (df["Cost"].between(thresholds["Cost"][0], thresholds["Cost"][1])) &
        (df["Solubility"].between(thresholds["Solubility"][0], thresholds["Solubility"][1])) &
        (df["Thermostability"].between(thresholds["Thermostability"][0], thresholds["Thermostability"][1])) &
        (df["Num_Alpha_Helices"].between(thresholds["Num_Alpha_Helices"][0], thresholds["Num_Alpha_Helices"][1])) &
        (df["Entropy"].between(thresholds["Entropy"][0], thresholds["Entropy"][1])) &
        (df["DNA_Complexity"].between(thresholds["DNA_Complexity"][0], thresholds["DNA_Complexity"][1]))
    ]
    return good_sequences

# Main logic to adjust thresholds to meet the required number of sequences
def adjust_thresholds(df, user_weightings, num_sequences):
    adjusted_thresholds = user_weightings.copy()
    increment = 1  # Adjust this value to control how aggressively the thresholds are lowered (Default = 1)

    while True:
        thresholds = {
            "Molecular_Weight": determine_thresholds_percentile(df, "Molecular_Weight", percentile=adjusted_thresholds["Molecular_Weight"]),
            "ASA": determine_thresholds_percentile(df, "ASA", percentile=adjusted_thresholds["ASA"]),
            "Isoelectric_Point": determine_thresholds_percentile(df, "Isoelectric_Point", percentile=adjusted_thresholds["Isoelectric_Point"]),
            "Cost": determine_thresholds_percentile(df, "Cost", percentile=adjusted_thresholds["Cost"]),
            "Solubility": determine_thresholds_percentile(df, "Solubility", percentile=100 - adjusted_thresholds["Solubility"], reverse=True),
            "Thermostability": determine_thresholds_percentile(df, "Thermostability", percentile=100 - adjusted_thresholds["Thermostability"], reverse=True),
            "Num_Alpha_Helices": determine_thresholds_percentile(df, "Num_Alpha_Helices", percentile=100 - adjusted_thresholds["Num_Alpha_Helices"], reverse=True),
            "Entropy": determine_thresholds_percentile(df, "Entropy", percentile=100 - adjusted_thresholds["Entropy"], reverse=True),
            "DNA_Complexity": determine_thresholds_percentile(df, "DNA_Complexity", percentile=100 - adjusted_thresholds["DNA_Complexity"], reverse=True)
        }

        selected_sequences = select_sequences(df, thresholds, num_sequences)

        if len(selected_sequences) >= num_sequences:
            # Trim the excess sequences
            sorted_sequences = selected_sequences.sort_values(
                by=["Molecular_Weight", "ASA", "Isoelectric_Point", "Cost", "Solubility", "Thermostability", "Num_Alpha_Helices", "Entropy", "DNA_Complexity"],
                ascending=[True, True, True, True, False, False, False, False, False]
                # Sorting in ascending order for the first four parameters to get the lowest values at the top.
                # Sorting in ascending order for the last five parameters to get the highest values at the top.
            )
            trimmed_sequences = sorted_sequences.head(num_sequences)
            return trimmed_sequences, thresholds

        for key in adjusted_thresholds.keys():
            if adjusted_thresholds[key] - increment > 0:
                adjusted_thresholds[key] -= increment

# User input prompt for number of sequences to select
num_sequences = int(input("How many sequences do you need? "))

# Prompt for solubility measurements
print("Run solubility measurements? (predict.py must be in present directory)")
print("Example command line execution: python predict.py --FASTA_PATH /home/s_gus/progs/D.fasta --OUTPUT_PATH /home/s_gus/progs/test_preds.csv --MODEL_TYPE ESM12 --PREDICTION_TYPE S")
run_solubility = input("Run solubility measurements? (Y/N): ").strip().upper()

if run_solubility == 'Y':
# If a prompt is wanted:
#    fasta_path = input("Enter the path for the FASTA file (e.g., /home/s_gus/progs/D.fasta): ").strip()
#    output_path = input("Enter the path for the output CSV file (e.g., /home/s_gus/progs/): ").strip()
    fasta_path = '/home/s_gus/progs/D.fasta'
    output_path = '/home/s_gus/progs/'
    thermo_path = '/home/s_gus/progs/long_sequence_predictions.tsv'
    if os.path.isfile('predict.py'):
        try:
            command = [
                'python', 'predict.py', 
                '--FASTA_PATH', fasta_path, 
                '--OUTPUT_PATH', output_path, 
                '--MODEL_TYPE', 'ESM12', 
                '--PREDICTION_TYPE', 'S'
            ]
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while running predict.py: {e}")
    else:
        print("predict.py not found in the present directory.")
else:
    print("Skipping solubility measurements.")
#    fasta_path = input("Enter the path for the FASTA file (e.g., /home/s_gus/progs/D.fasta): ").strip()
#    output_path = input("Enter the path for the output CSV file (e.g., /home/s_gus/progs/): ").strip()
    fasta_path = '/home/s_gus/progs/D.fasta'
    output_path = '/home/s_gus/progs/'    
    thermo_path = "/home/s_gus/progs/long_sequence_predictions.tsv"

# Solubility calculation from NETSOLP output
solubility_path = os.path.join(output_path, "test_preds.csv")
df = calculate_properties(fasta_path, solubility_path, thermo_path)
#last_thermostability_threshold = df['Thermostability']
#print(f"Last thermostability threshold: {last_thermostability_threshold}")
user_weightings = ask_user_for_weights()

# Adjust thresholds to meet the required number of sequences
selected_sequences, final_thresholds = adjust_thresholds(df, user_weightings, num_sequences)

# Write selected sequences to a FASTA file
fasta_output_path = os.path.join(output_path, "selected_sequences.fasta")
selected_sequences_list = [(seq_id, sequence) for seq_id, sequence in read_fasta(fasta_path) if seq_id in selected_sequences["ID"].values]
write_fasta(selected_sequences_list, fasta_output_path)

print(f"Selected {len(selected_sequences)} sequences written to {fasta_output_path}.")