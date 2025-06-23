#!/usr/bin/env python3

import argparse
import subprocess
import os
import shutil
from Bio import SeqIO
import pandas as pd
import re

def run_miranda(mirna_db, utr_chunk_fa, miranda_bin):
    """Run miRanda on a UTR chunk and return its stdout lines."""
    result = subprocess.run(
        [miranda_bin, mirna_db, utr_chunk_fa],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"miRanda failed on {utr_chunk_fa}:\n{result.stderr}")
    return result.stdout.splitlines()

def parse_miranda_output(lines):
    """Extract hits from verbose miRanda output."""
    hits = []
    hit_pattern = re.compile(
        r"^>(\S+)\s+(\S+)\s+([\d.]+)\s+(-?[\d.]+)\s+\d+\s+\d+\s+(\d+)\s+(\d+)"
    )
    for line in lines:
        match = hit_pattern.match(line)
        if match:
            mirna, query, score, energy, start, end = match.groups()
            hits.append({
                "query": query,
                "mirna": mirna,
                "score": float(score),
                "energy": float(energy),
                "start": int(start),
                "end": int(end)
            })
    return hits

def filter_overlapping_hits(hits):
    """Keep only best (lowest-energy) hit per overlapping region."""
    sorted_hits = sorted(hits, key=lambda x: (x['start'], x['end'], x['energy']))
    filtered = []
    current_group = []

    for hit in sorted_hits:
        if not current_group:
            current_group.append(hit)
            continue

        last = current_group[-1]
        if hit['start'] <= last['end']:
            current_group.append(hit)
        else:
            best = min(current_group, key=lambda x: x['energy'])
            filtered.append(best)
            current_group = [hit]

    if current_group:
        best = min(current_group, key=lambda x: x['energy'])
        filtered.append(best)

    return filtered

def write_output_table(hits, out_tsv):
    df = pd.DataFrame(hits)
    df = df[["query", "mirna", "score", "energy", "start", "end"]]
    df.to_csv(out_tsv, sep="\t", index=False, header=True)

def clean_temp_dir(temp_dir, prefix="utr_chunk_"):
    """Remove chunk FASTA files and other temp results."""
    for fname in os.listdir(temp_dir):
        if fname.startswith(prefix):
            try:
                os.remove(os.path.join(temp_dir, fname))
            except Exception as e:
                print(f"Warning: Failed to delete {fname}: {e}")

def main():
    parser = argparse.ArgumentParser(description="miRanda Python Wrapper (Verbose Output Compatible)")
    parser.add_argument("mirna_db", help="Path to miRNA database (FASTA)")
    parser.add_argument("utr_fasta", help="Path to UTR FASTA file")
    parser.add_argument("raw_output", help="Path to save raw miRanda output (combined)")
    parser.add_argument("--miranda_bin", required=True, help="Path to miRanda binary")
    parser.add_argument("--parsed_out", required=True, help="Path to save parsed output TSV")
    parser.add_argument("--tmpdir", required=True, help="Temporary directory for chunking")
    parser.add_argument("--chunk_size", type=int, default=500, help="Number of sequences per chunk")

    args = parser.parse_args()

    all_hits = []
    all_raw_lines = []
    records = list(SeqIO.parse(args.utr_fasta, "fasta"))

    if not os.path.exists(args.tmpdir):
        os.makedirs(args.tmpdir)

    # Process in chunks
    for i in range(0, len(records), args.chunk_size):
        chunk = records[i:i + args.chunk_size]
        chunk_path = os.path.join(args.tmpdir, f"utr_chunk_{i}.fa")
        SeqIO.write(chunk, chunk_path, "fasta")

        try:
            lines = run_miranda(args.mirna_db, chunk_path, args.miranda_bin)
            hits = parse_miranda_output(lines)
            all_raw_lines.extend(lines)
            all_hits.extend(hits)
        except Exception as e:
            print(f"Error in chunk {chunk_path}: {e}")

    # Save raw output (for debugging/record)
    with open(args.raw_output, "w") as raw_f:
        raw_f.write("\n".join(all_raw_lines))

    # Filter best hits
    best_hits = filter_overlapping_hits(all_hits)

    # Write final filtered hits
    write_output_table(best_hits, args.parsed_out)

    # Clean up
    clean_temp_dir(args.tmpdir)

if __name__ == "__main__":
    main()

