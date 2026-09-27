# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "biopython",
#     "numpy",
# ]
# ///

from Bio.Seq import Seq
from Bio.Data import CodonTable
import numpy as np
import sys
import os

# A simplified E. coli codon bias table (most frequent codons for each amino acid)
E_COLI_CODONS = {
    'A': 'GCG', 'R': 'CGC', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC',
    'Q': 'CAG', 'E': 'GAA', 'G': 'GGC', 'H': 'CAT', 'I': 'ATT',
    'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCG',
    'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTG',
    '*': 'TAA' # Stop codon
}

def reverse_translate(protein_seq):
    """Converts Amino Acid sequence back into optimal DNA sequence for E. coli"""
    dna = []
    for aa in protein_seq:
        if aa in E_COLI_CODONS:
            dna.append(E_COLI_CODONS[aa])
        else:
            raise ValueError(f"Unknown amino acid: {aa}")
    return "".join(dna)

def generate_orca_input(fasta_path, output_dir):
    print(f"Reading Protein FASTA: {fasta_path}")
    with open(fasta_path, "r") as f:
        lines = f.readlines()
        protein_seq = "".join([l.strip() for l in lines if not l.startswith(">")])
    
    print(f"Loaded TriFusion Protein ({len(protein_seq)} Amino Acids)")
    
    # 1. Reverse Translation to DNA
    print("Reverse-Translating to optimized E. coli DNA sequence...")
    gene_dna = reverse_translate(protein_seq)
    
    # 2. Build the Synthetic Plasmid
    # We add a generic T7 Promoter (to start transcription) and a T7 Terminator
    T7_PROMOTER = "TAATACGACTCACTATAGGGGGAATTGTGAGCGGATAACAATTCCCCT"
    T7_TERMINATOR = "TAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG"
    
    synthetic_plasmid = T7_PROMOTER + gene_dna + E_COLI_CODONS['*'] + T7_TERMINATOR
    
    print(f"Constructed Synthetic DNA Plasmid ({len(synthetic_plasmid)} base pairs)")
    
    # 3. Format for ORCA (Sequence-based modeling of 3D genome architecture)
    # ORCA expects genomic sequences in specific FASTA/bed formats to predict the 
    # hierarchical chromatin folding (TADs, loops) via its CNN architecture.
    os.makedirs(output_dir, exist_ok=True)
    
    orca_input_path = os.path.join(output_dir, "teaflon_orca_input.fasta")
    with open(orca_input_path, "w") as f:
        f.write(">chr_synthetic_teaflon_plasmid\n")
        f.write(synthetic_plasmid + "\n")
        
    print(f"\nSaved ORCA DNA Input to: {orca_input_path}")
    print("\n--- ORCA Deep Learning Validation Step ---")
    print("The synthetic DNA sequence is now ready to be fed into the ORCA model.")
    print("ORCA will use its Convolutional Neural Networks to predict the 3D chromatin folding.")
    print("Validation Goal: Ensure the T7 Promoter physically loops to the Gene Start Codon")
    print("and that the synthetic sequence does not form a tightly closed heterochromatin knot")
    print("that would prevent E. coli from transcribing our TeaFlon enzyme!")

if __name__ == "__main__":
    fasta_file = "teaflon_trifusion.fasta"
    if not os.path.exists(fasta_file):
        print(f"Error: {fasta_file} not found.")
        sys.exit(1)
        
    generate_orca_input(fasta_file, "orca_validation")
