#!/usr/local/bin/python3

# Read in the genome of the virus
genomeFile = open('sars-cov-2_complete-genome-2019-12', 'r')
sars_cov_2_g0 = genomeFile.read()

# Determine amount of data in the virus
print(f"SARS-COV-2 genome is {len(sars_cov_2_g0)} characters long, meaning it has {len(sars_cov_2_g0) / 1000}kB of data")

# Now translate into protein sequence

# store the codons to amino acid table
table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

# check that the DNA sequence is divisible by 3
print(sars_cov_2_g0)
if len(sars_cov_2_g0) % 3 != 0:
    print("[!] DNA Sequence not divisible by 3")
    exit(1)
