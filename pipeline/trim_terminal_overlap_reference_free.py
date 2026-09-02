#!/usr/bin/env python3

import argparse
import subprocess
import tempfile
from pathlib import Path

def read_fasta(path):
    name = None
    seqs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
            else:
                seqs.append(line)
    return name or "contig", "".join(seqs)

def write_fasta(path, name, seq):
    with open(path, "w") as out:
        out.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")

def parse_paf(paf):
    hits = []
    with open(paf) as f:
        for line in f:
            p = line.rstrip().split("\t")
            if len(p) < 12:
                continue

            qname = p[0]
            qlen = int(p[1])
            qs = int(p[2])
            qe = int(p[3])
            strand = p[4]
            tname = p[5]
            tlen = int(p[6])
            ts = int(p[7])
            te = int(p[8])
            matches = int(p[9])
            block = int(p[10])
            mapq = int(p[11])
            identity = matches / block if block else 0

            hits.append({
                "qname": qname,
                "qlen": qlen,
                "qs": qs,
                "qe": qe,
                "strand": strand,
                "tname": tname,
                "tlen": tlen,
                "ts": ts,
                "te": te,
                "matches": matches,
                "block": block,
                "mapq": mapq,
                "identity": identity,
            })
    return hits

def main():
    parser = argparse.ArgumentParser(
        description="Reference-free terminal redundancy/circular overlap trimming for a single-contig FASTA."
    )
    parser.add_argument("input_fasta")
    parser.add_argument("output_fasta")
    parser.add_argument("--window", type=int, default=700000,
                        help="Number of bases from each end to compare")
    parser.add_argument("--min-overlap", type=int, default=10000,
                        help="Minimum terminal overlap length to trim")
    parser.add_argument("--min-identity", type=float, default=0.99,
                        help="Minimum alignment identity")
    parser.add_argument("--terminal-tolerance", type=int, default=2000,
                        help="Allowed distance from exact contig end/start")
    parser.add_argument("--minimap2", default="minimap2")
    args = parser.parse_args()

    name, seq = read_fasta(args.input_fasta)
    L = len(seq)

    if L < 2 * args.min_overlap:
        raise SystemExit("ERROR: sequence too short for terminal-overlap trimming")

    w = min(args.window, L // 2)

    prefix = seq[:w]
    suffix = seq[-w:]

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        prefix_fa = tmp / "prefix.fa"
        suffix_fa = tmp / "suffix.fa"
        paf = tmp / "suffix_vs_prefix.paf"

        write_fasta(prefix_fa, "prefix", prefix)
        write_fasta(suffix_fa, "suffix", suffix)

        cmd = [
            args.minimap2,
            "-x", "asm5",
            str(prefix_fa),
            str(suffix_fa)
        ]

        with open(paf, "w") as out:
            subprocess.run(cmd, stdout=out, check=True)

        hits = parse_paf(paf)

    candidates = []

    for h in hits:
        aln_len = h["qe"] - h["qs"]

        if h["identity"] < args.min_identity:
            continue

        if aln_len < args.min_overlap:
            continue

        if h["strand"] != "+":
            continue

        # Case 1:
        # end of suffix aligns to start of prefix.
        # This means the assembly end repeats the assembly beginning.
        suffix_touches_end = abs(h["qe"] - h["qlen"]) <= args.terminal_tolerance
        prefix_touches_start = h["ts"] <= args.terminal_tolerance

        if suffix_touches_end and prefix_touches_start:
            candidates.append({
                "action": "trim_end",
                "trim_len": aln_len,
                "identity": h["identity"],
                "block": h["block"],
                "details": h
            })

        # Case 2:
        # start of suffix aligns to end of prefix.
        # Less common here, but included for completeness.
        suffix_touches_start = h["qs"] <= args.terminal_tolerance
        prefix_touches_end = abs(h["te"] - h["tlen"]) <= args.terminal_tolerance

        if suffix_touches_start and prefix_touches_end:
            candidates.append({
                "action": "trim_start",
                "trim_len": aln_len,
                "identity": h["identity"],
                "block": h["block"],
                "details": h
            })

    if not candidates:
        print("input_length:", L)
        print("No terminal duplicated overlap found. Writing unchanged sequence.")
        write_fasta(args.output_fasta, name + "_terminal_checked", seq)
        print("output_length:", L)
        print("saved:", args.output_fasta)
        return

    # Prefer the longest high-identity terminal overlap
    best = sorted(candidates, key=lambda x: (x["trim_len"], x["identity"]), reverse=True)[0]

    if best["action"] == "trim_end":
        trimmed = seq[:-best["trim_len"]]
    elif best["action"] == "trim_start":
        trimmed = seq[best["trim_len"]:]
    else:
        raise SystemExit("ERROR: unknown trimming action")

    print("input_length:", L)
    print("action:", best["action"])
    print("terminal_overlap_trimmed:", best["trim_len"])
    print("overlap_identity_pct:", best["identity"] * 100)
    print("output_length:", len(trimmed))
    print("alignment_details:", best["details"])

    write_fasta(args.output_fasta, name + "_terminal_overlap_trimmed", trimmed)
    print("saved:", args.output_fasta)

if __name__ == "__main__":
    main()
