#!/usr/bin/env python3
import argparse
import pandas as pd
from primer3 import bindings
parser = argparse.ArgumentParser(description="Design primers using primer3 for each gene in a FASTA file.")
parser.add_argument("input_fasta", help="Input FASTA file")
parser.add_argument("output_csv", help="Output CSV file")
parser.add_argument("--num_return", type=int, default=1, help="Number of primer pairs to return (default: 1)")
args = parser.parse_args()

def parse_fasta(fasta_path):
    records = []
    with open(fasta_path) as fh:
        seq_id, seq_parts = None, []
        for line in fh:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if seq_id:
                    records.append((seq_id, ''.join(seq_parts)))
                seq_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if seq_id:
            records.append((seq_id, ''.join(seq_parts)))
    return records

def design_primers(seq_id, seq, num_return):
    seq_args = {'SEQUENCE_ID': seq_id, 'SEQUENCE_TEMPLATE': seq}
    global_args = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_NUM_RETURN': num_return
    }
    results = []
    res = bindings.designPrimers(seq_args, global_args)
    for i in range(num_return):
        results.append({
            "gene_id": seq_id,
            "left_primer": res.get(f"PRIMER_LEFT_{i}_SEQUENCE", ""),
            "right_primer": res.get(f"PRIMER_RIGHT_{i}_SEQUENCE", ""),
            "left_tm": res.get(f"PRIMER_LEFT_{i}_TM", ""),
            "right_tm": res.get(f"PRIMER_RIGHT_{i}_TM", ""),
            "product_size": res.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", "")
        })
    return results

all_records = parse_fasta(args.input_fasta)
all_primers = []
for seq_id, seq in all_records:
    all_primers.extend(design_primers(seq_id, seq, args.num_return))

pd.DataFrame(all_primers).to_csv(args.output_csv, index=False)
