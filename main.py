#!/usr/local/bin/python3

import requests
import re

# Functions

def translate(rna):

    # store the codons to amino acid table
    aa_table = {
        'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
        'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
        'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
        'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
        'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
        'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
        'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
        'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
        'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
        'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W',
    }

    # check that the RNA sequence is divisible by 3
    if len(rna) % 3 != 0:
        print("[!] RNA Sequence not divisible by 3")
        # exit(1)

    # split the rna sequence string into codons and store them into an array
    sars_cov_2_codons_arr = []
    counter = 0
    while counter < len(rna):
        sars_cov_2_codons_arr.append(rna[counter:counter+3])
        counter += 3

    # now Translate codon by codon and return the amino acid sequence string
    sars_cov_2_aa_raw = []
    for codon in sars_cov_2_codons_arr:
        sars_cov_2_aa_raw.append(aa_table[codon])

    sars_cov_2_aa_raw= "".join(sars_cov_2_aa_raw)
    return sars_cov_2_aa_raw


# program start
# Read in the genome of the virus
genomeFile = open('MN908947.3_sars-cov-2.fasta', 'r')
sars_cov_2_g0_raw = genomeFile.read().rstrip().split("\n")

sars_cov_2_g0 = {
    "description": "",
    "genome": "",
    "rna": ""
}

for line in sars_cov_2_g0_raw:
    if line[0] == ">":
        sars_cov_2_g0["description"] = line
    else:
        sars_cov_2_g0["genome"] += line

# print(sars_cov_2_g0) #debug

# Determine amount of data in the virus
print(f"SARS-COV-2 genome is {len(sars_cov_2_g0['genome'])} characters long, meaning it has {len(sars_cov_2_g0['genome']) / 1000}kB of data")
print("----------------------------")
## Now translate into rna

sars_cov_2_g0['rna'] = sars_cov_2_g0['genome'].replace("T", "U")

# store parts of the genomic rna in variables
sars_cov_2_5utr = sars_cov_2_g0['rna'][0:265]
sars_cov_2_transl = sars_cov_2_g0['rna'][265:29674]
sars_cov_2_3utr = sars_cov_2_g0['rna'][29674:29903]

NCBI_page = requests.get('https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=MN908947.3&db=nuccore&report=genbank&conwithfeat=on&withparts=on&hide-cdd=on&retmode=text&withmarkup=on&tool=portal&log$=seqview&maxdownloadsize=1000000')

split_data_ncbi = NCBI_page.text.split("\n")
sars_cov_2_genes = []
gene = {}

# go thru the scraped page line by line and put the gene data in order
for line in split_data_ncbi:
    # get location
    geneloc_regex = r"\d+\.\.\d+"
    location_matches = re.findall(geneloc_regex, line, re.MULTILINE)
    if "CDS" in line and len(location_matches) > 0:
        gene_location = []
        for match in location_matches:
            location_split = match.split("..")
            location_parsed = {
             "from": int(location_split[0]),
             "to": int(location_split[1])
            }
            gene_location.append(location_parsed)

        gene["location"] = gene_location


    # get gene name
    if "/gene=" in line:
        gene["gene_name"] = line.split("\"")[1]
    # get protein name
    if "/product=" in line:
        gene["protein_name"] = line.split("\"")[1]
    # get protein_id
    if "/protein_id=" in line:
        gene["protein_id"] = line.split("\"")[1]
        # since protein_id is the last data point for each gene we
        # now get the rna part
        rna_seq = ""
        for coordinates in gene["location"]:
            rna_from = coordinates["from"] - 1
            rna_to = coordinates["to"]
            rna_seq += sars_cov_2_g0['rna'][rna_from:rna_to]

        gene["rna"] = rna_seq
        # add the gene object to the genes array and empty it out
        # to free space for the next gene
        sars_cov_2_genes.append(gene)
        gene = {}


for gene in sars_cov_2_genes:
    gene['aa'] = translate(gene['rna'])
    print("-----")
    print(gene)
    print("-----")


exit(0)
