# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "torch",
#     "transformers",
# ]
# ///

import torch
from transformers import AutoTokenizer, EsmForProteinFolding
import sys

def fold_sequence(sequence, output_pdb_path):
    print("Loading ESMFold tokenizer...")
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    
    print("Loading ESMFold model (This will download ~11GB of weights on the first run)...")
    # Use MPS (Apple Silicon GPU) if available, otherwise CPU
    device = "mps" if torch.backends.mps.is_available() else "cpu"
    print(f"Using device: {device}")
    
    # Load model with low_cpu_mem_usage to help with RAM during load
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    model = model.to(device)
    
    print(f"Tokenizing sequence ({len(sequence)} amino acids)...")
    inputs = tokenizer([sequence], return_tensors="pt", add_special_tokens=False)
    inputs = {k: v.to(device) for k, v in inputs.items()}
    
    print("Running structure prediction... (This may take a few minutes)")
    with torch.no_grad():
        outputs = model(**inputs)
    
    print("Prediction complete! Formatting to PDB...")
    # Convert tensor outputs to PDB string format
    pdb_string = model.output_to_pdb(outputs)[0]
    
    with open(output_pdb_path, "w") as f:
        f.write(pdb_string)
        
    print(f"Successfully saved folded structure to {output_pdb_path}")

if __name__ == "__main__":
    # Read our TriFusion fasta
    with open("teaflon_trifusion.fasta", "r") as f:
        lines = f.readlines()
        # Join lines ignoring the header
        seq = "".join([l.strip() for l in lines if not l.startswith(">")])
        
    fold_sequence(seq, "teaflon_trifusion_local.pdb")
