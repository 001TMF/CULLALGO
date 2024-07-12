# CULLALGO
'CULLALGO is a script that uses multiple parameters to 'cull' a large dataset of .fasta protein sequences to a manageable amount with desired traits in a .fasta output. All measurements are based on existing literature.'

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

'NETSOLP https://academic.oup.com/bioinformatics/article/38/4/941/6444984'

'TemStaPro https://github.com/ievapudz/TemStaPro'

# Set-up

1. Install and correctly set-up NETSOLP and TemStaPro. It is advised to first test both programs independently to make sure they work. As of CULLALGO1.0, NETSOLP Sollubility
measurements CAN be run natively in the script after set-up. However, TemStaPro MUST BE FIRST RUN SEPERATELY, and the output file correctly placed in the desired directory. It
cannot be run natively in the script.
2. 
