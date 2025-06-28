#!/usr/bin/env python3

import argparse, os, subprocess, re
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

def split_into_chunks(seq_record, chunk_size, step_size, tmpdir, job_id):
    seq = str(seq_record.seq)
    chunks = []
    for start in range(0, len(seq), step_size):
        end = min(start + chunk_size, len(seq))
        chunk_seq = seq[start:end]
        chunk_name = f"{job_id}_chunk_{start+1}_{end}"
        chunk_path = os.path.join(tmpdir, f"{chunk_name}.fa")
        record = SeqRecord(Seq(chunk_seq), id=chunk_name, description="")
        SeqIO.write(record, chunk_path, "fasta")
        chunks.append((chunk_path, start + 1))
        if end == len(seq):
            break
    return chunks

def run_miranda(mirna_db, chunk_path, miranda_bin, offset, min_energy, min_score, require_seed_start):
    try:
        result = subprocess.run(
            [miranda_bin, mirna_db, chunk_path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        lines = result.stdout.splitlines()
        return parse_hits(lines, offset, min_energy, min_score, require_seed_start), result.stdout
    except Exception as e:
        return [], f"# ERROR: {chunk_path}: {e}"

def parse_hits(lines, offset, min_energy, min_score, require_seed_start):
    hits = []
    hit_pattern = re.compile(
        r"^>(\S+)\s+(\S+)\s+([\d.]+)\s+(-?[\d.]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)"
    )
    for line in lines:
        match = hit_pattern.match(line)
        if match:
            mirna, query, score, energy, qstart, qend, tstart, tend = match.groups()
            score, energy, qstart = float(score), float(energy), int(qstart)

            if score < min_score or energy > min_energy:
                continue
            if require_seed_start and qstart > 2:
                continue

            hits.append({
                "query": query,
                "mirna": mirna,
                "score": score,
                "energy": energy,
                "start": offset + int(tstart) - 1,
                "end": offset + int(tend) - 1
            })
    return hits

def filter_hits_by_overlap(hits, min_overlap_pct):
    sorted_hits = sorted(hits, key=lambda x: (x["start"], x["end"], x["energy"]))
    filtered = []
    current_group = []

    def overlap_pct(a, b):
        ovl = max(0, min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1)
        union = max(a["end"], b["end"]) - min(a["start"], b["start"]) + 1
        return ovl / union

    for hit in sorted_hits:
        if not current_group:
            current_group.append(hit)
        elif overlap_pct(current_group[-1], hit) >= min_overlap_pct:
            current_group.append(hit)
        else:
            best = min(current_group, key=lambda x: x["energy"])
            filtered.append(best)
            current_group = [hit]
    if current_group:
        best = min(current_group, key=lambda x: x["energy"])
        filtered.append(best)
    return filtered

def save_hits_tsv(hits, out_path):
    df = pd.DataFrame(hits)[["query", "mirna", "score", "energy", "start", "end"]]
    df.to_csv(out_path, sep="\t", index=False)

def cleanup_chunks(tmpdir, job_id):
    for f in os.listdir(tmpdir):
        if f.startswith(job_id) and f.endswith(".fa"):
            os.remove(os.path.join(tmpdir, f))

def main():
    parser = argparse.ArgumentParser(description="miRanda wrapper with score, energy, seed, and overlap filters.")
    parser.add_argument("mirna_db", help="miRNA database FASTA")
    parser.add_argument("utr_fasta", help="Single UTR target FASTA")
    parser.add_argument("raw_output", help="Raw combined miRanda output path")
    parser.add_argument("--parsed_out", required=True, help="TSV output for filtered hits")
    parser.add_argument("--miranda_bin", required=True, help="Path to miRanda binary")
    parser.add_argument("--tmpdir", required=True, help="Temporary directory for chunk files")
    parser.add_argument("--job_id", required=True, help="Job ID for chunk file naming")
    parser.add_argument("--chunk_size", type=int, default=500, help="Window size for sequence chunking")
    parser.add_argument("--step_size", type=int, default=400, help="Step size between chunks")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads")
    parser.add_argument("--min_energy", type=float, default=-25.0, help="Minimum allowed energy")
    parser.add_argument("--min_score", type=float, default=155.0, help="Minimum miRanda alignment score")
    parser.add_argument("--min_overlap_pct", type=float, default=0.8, help="Minimum overlap % to merge hits")
    parser.add_argument("--require_seed_start", action='store_true', default=True,
                        help="Require miRNA alignment to start at 5′ end (default: True)")

    args = parser.parse_args()
    os.makedirs(args.tmpdir, exist_ok=True)

    records = list(SeqIO.parse(args.utr_fasta, "fasta"))
    if len(records) != 1:
        raise ValueError("UTR FASTA must contain exactly one sequence.")

    chunks = split_into_chunks(records[0], args.chunk_size, args.step_size, args.tmpdir, args.job_id)

    all_hits, all_raw = [], []

    with ProcessPoolExecutor(max_workers=args.threads) as pool:
        futures = {
            pool.submit(run_miranda, args.mirna_db, path, args.miranda_bin, offset,
                        args.min_energy, args.min_score, args.require_seed_start): path
            for path, offset in chunks
        }

        for future in as_completed(futures):
            hits, raw = future.result()
            all_hits.extend(hits)
            all_raw.append(raw)

    with open(args.raw_output, "w") as raw_out:
        raw_out.write("\n".join(all_raw))

    filtered = filter_hits_by_overlap(all_hits, args.min_overlap_pct)
    save_hits_tsv(filtered, args.parsed_out)
    cleanup_chunks(args.tmpdir, args.job_id)

if __name__ == "__main__":
    main()
