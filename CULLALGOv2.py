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
import yaml
import argparse
from tqdm import tqdm
import time

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
# Adjust window size if needed (default = 10)----------------------------------------------------------------------
# Adjust threshold for narrower picking (default = 1.0)------------------------------------------------------------
def calculate_propensity(sequence, window_size=10, threshold=1.0):
    propensities = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i + window_size]
        avg_propensity = sum(alpha_helix_propensity.get(aa, 0) for aa in window) / window_size
        propensities.append(avg_propensity)
    return propensities


# Iterates propensities to determine helical regions
# Window size and threshold should match from calculate_propensity-------------------------------------------------
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


def count_sequences_in_fasta(fasta_path):
    """Counts the number of sequences in a given FASTA file."""
    count = 0
    with open(fasta_path, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            count += 1
    return count


# Function to parse command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Process command line arguments.')
    parser.add_argument('--config', type=str, default='config.yaml', help='Path to the configuration file')
    return parser.parse_args()


# Function to read configuration from YAML file
def load_config(file_path):
    with open(file_path, 'r') as file:
        return yaml.safe_load(file)


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


def read_thermostability(thermo_path):
    df = pd.read_csv(thermo_path, sep='\t')
    df.columns = [f'Col{i + 1}' for i in range(len(df.columns))]
    thermostability_column = df.columns[16]
    df['Thermostability'] = df[thermostability_column].apply(extract_last_number)
    return df[['Col1', 'Thermostability']].rename(columns={'Col1': 'ID'})


def extract_last_number(value):
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
    entropy = -sum((count / length) * math.log2(count / length) for count in counts.values())
    return entropy


# Calculates all properties in format
def calculate_properties(fasta_path, solubility_path, thermo_path):
    sequences = read_fasta(fasta_path)
    solubility_data = read_solubility(solubility_path)
    thermostability_data = read_thermostability(thermo_path)

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
            thermostability = thermostability_row['Thermostability'].values[
                0] if not thermostability_row.empty else None

        data.append((seq_id, weight, asa, isoelectric_point, cost, solubility, thermostability, num_alpha_helices,
                     entropy, dna_complexity))

    df = pd.DataFrame(data, columns=['ID', 'Molecular_Weight', 'ASA', 'Isoelectric_Point', 'Cost', 'Solubility',
                                     'Thermostability', 'Num_Alpha_Helices', 'Entropy', 'DNA_Complexity'])
    return df


# Threshold percentages for manual tweaking of starting thresholds
# Percentile dictates the pool of sequences the script picks from. If percentile=100, it has access to all sequences.
# I.e., the minimum is the 0th percentile.
# If e.g., percentile=80, the minimum is the 20th percentile, so the bottom 20% is not accessible.
# It is advised to leave this unchanged for best results (default percentile=100)---------------------------------------------------------------------
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
# The percentiles associated with weight integers [1,2,3] are default as:
# 1 - 0th percentile. All sequences included.
# 2 - 50th percentile. 50% of the worst sequences culled.
# 3 - 80th percentile. 80% of the worst sequences culled.
# weightings[parameter] can be changed according to needs-------------------------------------------------------------------------
def ask_user_for_weights(config):
    parameters = ["Molecular_Weight", "ASA", "Isoelectric_Point", "Cost", "Solubility", "Thermostability",
                  "Num_Alpha_Helices", "Entropy", "DNA_Complexity"]
    weightings = {}
    for parameter in parameters:
        # Read the weighting from the config file
        weight = config['weights'].get(parameter, 2)  # default to 2 if not specified
        if weight == 1:
            weightings[parameter] = 0
        elif weight == 2:
            weightings[parameter] = 50
        else:
            weightings[parameter] = 100
    return weightings


def write_fasta(sequences, output_path):
    fasta_records = []
    for seq_id, sequence in sequences:
        record = SeqRecord(Seq(sequence), id=seq_id, description="")
        fasta_records.append(record)
    SeqIO.write(fasta_records, output_path, "fasta")


# Select sequences from user thresholds
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
def adjust_thresholds(df, user_weightings, num_sequences, max_iterations=1000):
    adjusted_thresholds = user_weightings.copy()
    increment = 1  # Adjust this value to control how aggressively the thresholds are lowered (Default = 1)

    iteration = 0

    while iteration < max_iterations:
        thresholds = {
            "Molecular_Weight": determine_thresholds_percentile(df, "Molecular_Weight", percentile=adjusted_thresholds["Molecular_Weight"]),
            "ASA": determine_thresholds_percentile(df, "ASA", percentile=adjusted_thresholds["ASA"]),
            "Isoelectric_Point": determine_thresholds_percentile(df, "Isoelectric_Point", percentile=adjusted_thresholds["Isoelectric_Point"]),
            "Cost": determine_thresholds_percentile(df, "Cost", percentile=adjusted_thresholds["Cost"]),
            "Solubility": determine_thresholds_percentile(df, "Solubility", percentile=100 - adjusted_thresholds["Solubility"], reverse=True),
            "Thermostability": determine_thresholds_percentile(df, "Thermostability", percentile=100 - adjusted_thresholds["Thermostability"], reverse=True),
            "Num_Alpha_Helices": determine_thresholds_percentile(df, "Num_Alpha_Helices", percentile=100 - adjusted_thresholds["Num_Alpha_Helices"], reverse=True),
            "Entropy": determine_thresholds_percentile(df, "Entropy", percentile=adjusted_thresholds["Entropy"]),  # Adjusted to filter for lowest entropy
            "DNA_Complexity": determine_thresholds_percentile(df, "DNA_Complexity", percentile=100 - adjusted_thresholds["DNA_Complexity"], reverse=True)
        }

        selected_sequences = select_sequences(df, thresholds, num_sequences)

        if len(selected_sequences) >= num_sequences:
            # Trim the excess sequences
            sorted_sequences = selected_sequences.sort_values(
                by=["Molecular_Weight", "ASA", "Isoelectric_Point", "Cost", "Solubility", "Thermostability", "Num_Alpha_Helices", "Entropy", "DNA_Complexity"],
                ascending=[True, True, True, True, False, False, False, True, False]
                # Sorting in ascending order for the first four parameters to get the lowest values at the top.
                # Sorting in ascending order for Entropy to get the lowest values at the top.
                # Sorting in ascending order for the last parameters to get the highest values at the top.
            )
            trimmed_sequences = sorted_sequences.head(num_sequences)
            return trimmed_sequences, thresholds

        for key in adjusted_thresholds.keys():
            if adjusted_thresholds[key] - increment > 0:
                adjusted_thresholds[key] -= increment

        iteration += 1

    raise ValueError("Unable to find the required number of sequences within the given constraints.")


def match_files(fasta_path, thermo_path):
    fasta_files = {os.path.splitext(f)[0]: f for f in os.listdir(fasta_path) if f.endswith('.fasta')}
    thermo_files = {os.path.splitext(f)[0].replace('-TemStaPro', ''): f for f in os.listdir(thermo_path) if f.endswith('.tsv')}
    return [(os.path.join(fasta_path, fasta_files[name]), os.path.join(thermo_path, thermo_files[name])) for name in thermo_files if name in fasta_files]

def main():
    args = parse_args()
    config = load_config(args.config)

    fasta_directory = config['paths']['fasta']
    output_path = config['paths']['output']
    thermo_directory = config['paths']['thermo']
    run_solubility = config.get('run_solubility', 'Y')
    culled_fasta = 'culled_fasta'

    matched_files = match_files(fasta_directory, thermo_directory)
    total_files = len(matched_files)

    pbar = tqdm(total=total_files, desc='Processing FASTA Files', unit='file')
    for fasta_file, thermo_file in matched_files:
        start_time = time.time()
        fasta_data = read_fasta(fasta_file)
        total_sequences = len(fasta_data)
        num_sequences = config.get("num_sequences", total_sequences)

        print(f"You are culling down to {num_sequences} sequence(s) out of a maximum {total_sequences}")
        if not 1 <= num_sequences <= total_sequences:
            raise ValueError(f"Number of sequences must be between 1 and {total_sequences}")

        print(f"{fasta_file} culling has begun")

        base_name = os.path.splitext(os.path.basename(fasta_file))[0]
        csv_output_name = f"{base_name}_solubility.csv"
        solubility_path = os.path.join(output_path, csv_output_name)
        if not os.path.exists(solubility_path):
            open(solubility_path, 'w').close()

        if run_solubility == 'Y':
            if os.path.isfile('netsolp/predict.py'):
                command = [
                    'python', 'predict.py',
                    '--FASTA_PATH', fasta_file,
                    '--OUTPUT_PATH', solubility_path,
                    '--MODEL_TYPE', 'ESM12',
                    '--PREDICTION_TYPE', 'S'
                ]
                subprocess.run(command, check=True)
            else:
                print("predict.py not found in the present directory.")
        else:
            if 'external_solubility_path' in config:
                solubility_path = config['external_solubility_path']
                print("Using external solubility data from:", solubility_path)
            else:
                raise ValueError("External solubility path must be provided in the config when solubility calculations are skipped.")
        print(solubility_path)
        df = calculate_properties(fasta_file, solubility_path, thermo_file)
        user_weightings = ask_user_for_weights(config)

        selected_sequences, final_thresholds = adjust_thresholds(df, user_weightings, num_sequences)

        fasta_output_path = os.path.join(output_path, culled_fasta, f"{base_name}-Culled.fasta")
        selected_sequences_list = [(seq_id, sequence) for seq_id, sequence in fasta_data if seq_id in selected_sequences["ID"].values]
        write_fasta(selected_sequences_list, fasta_output_path)

        print(f"Selected {len(selected_sequences)} sequences written to {fasta_output_path}.")
        #update the progress bar
        pbar.update(1)
    pbar.close()


if __name__ == "__main__":
    main()
