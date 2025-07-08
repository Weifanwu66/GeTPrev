#!/usr/bin/env python3
import argparse
from Bio import SeqIO
parser = argparse.ArgumentParser(description="Extract sequences from a FASTA file based on a gene ID list.")
parser.add_argument("input_fasta", help="Input full FASTA file")
parser.add_argument("gene_list", help="Text file with one gene ID per line")
parser.add_argument("output_fasta", help="Output FASTA file with extracted genes")
args = parser.parse_args()
with open(args.gene_list) as f:
 target_ids = set(line.strip() for line in f if line.strip())
count = 0
with open(args.output_fasta, "w") as out_f:
 for record in SeqIO.parse(args.input_fasta, "fasta"):
  if record.id in target_ids:
   SeqIO.write(record, out_f, "fasta")
   count += 1
print(f"Extraction complete: {count} sequence written to '{args.output_fasta}'")
