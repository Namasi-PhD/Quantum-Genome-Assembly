#!/usr/bin/env python3
"""
merge_reports.py — merge per-sample QUAST report columns into one wide tsv.

Each file in results_columns/ is a copy of a QUAST report.tsv:
    Assembly          <sample>_polished_assembly
    # contigs (...)   1
    ...

This script combines them into a single tsv with one "Assembly" header row
listing every sample, and one row per metric label -- matching the layout
of report-real-device.tsv.

Usage:
    python3 merge_reports.py results_columns/ final_report.tsv
"""
import sys
import os
import glob
import re


def natural_sample_key(path):
    """Sort by the numeric sample id in the filename (0, 1, 2, ... 5000)."""
    name = os.path.splitext(os.path.basename(path))[0]
    m = re.match(r"(\d+)", name)
    return (int(m.group(1)) if m else float("inf"), name)


def read_report(path):
    header = None
    rows = {}
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            label = parts[0]
            value = parts[1] if len(parts) > 1 else ""
            if label == "Assembly":
                header = value
            else:
                rows[label] = value
    return header, rows


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <results_columns_dir> <output.tsv>")
        sys.exit(1)

    results_dir, out_path = sys.argv[1], sys.argv[2]
    files = sorted(glob.glob(os.path.join(results_dir, "*.tsv")), key=natural_sample_key)

    if not files:
        print(f"No .tsv files found in {results_dir}")
        sys.exit(1)

    all_labels = []
    seen_labels = set()
    columns = []  # list of (header, rows_dict)
    skipped = []

    for path in files:
        header, rows = read_report(path)
        if header is None:
            skipped.append(path)
            continue
        columns.append((header, rows))
        for label in rows:
            if label not in seen_labels:
                seen_labels.add(label)
                all_labels.append(label)

    with open(out_path, "w") as out:
        out.write("Assembly\t" + "\t".join(h for h, _ in columns) + "\n")
        for label in all_labels:
            values = [rows.get(label, "") for _, rows in columns]
            out.write(label + "\t" + "\t".join(values) + "\n")

    print(f"Merged {len(columns)} sample(s) into {out_path}")
    if skipped:
        print(f"WARNING: skipped {len(skipped)} malformed file(s):")
        for p in skipped:
            print(f"  - {p}")


if __name__ == "__main__":
    main()
