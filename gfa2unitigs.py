#!/usr/bin/env python3
import sys, re
from collections import defaultdict, deque

# ------------------------------------
# Parse the GFA
# ------------------------------------

segments = {}             # segment_name -> sequence
adj = defaultdict(list)   # oriented_end -> list of (next_oriented_end, cigar)

def revcomp(seq):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]

with open(sys.argv[1]) as f:
    for line in f:
        parts = line.strip().split("\t")
        if not parts:
            continue

        if parts[0] == "S":
            _, name, seq, *_ = parts
            segments[name] = seq

        elif parts[0] == "L":
            _, a, a_or, b, b_or, cigar, *_ = parts

            from_end = (a, a_or)
            to_end   = (b, b_or)

            # Add oriented edge
            adj[from_end].append((to_end, cigar))

            # Also add the reverse direction
            comp_from = (a, "-" if a_or == "+" else "+")
            comp_to   = (b, "-" if b_or == "+" else "+")
            adj[comp_to].append((comp_from, cigar))

# ------------------------------------
# Build undirected connectivity to find components
# ------------------------------------

undirected = defaultdict(set)
for (seg, ori), neighbors in adj.items():
    for (nseg, nori), _ in neighbors:
        undirected[seg].add(nseg)
        undirected[nseg].add(seg)

# ------------------------------------
# Walk a component to find path
# ------------------------------------

def find_component_start(component):
    """
    For a linear component with one path, find
    the segment end that should be the start
    based on degree.
    """
    end_degrees = defaultdict(int)
    for seg in component:
        for ori in ["+", "-"]:
            end_degrees[(seg, ori)] = len(adj[(seg, ori)])

    # Terminal ends have degree == 1
    candidates = [end for end, d in end_degrees.items() if d == 1 and end[0] in component]
    if not candidates:
        # fallback: take any seg+ as start
        return (next(iter(component)), "+")
    return candidates[0]

def walk_path(start_end, component):
    """
    Walk the oriented graph from start_end
    for a given component.
    Return list of (segment, orientation) and list of cigars.
    """
    visited = set()
    path = []
    cigars = []

    prev = None
    cur = start_end

    while True:
        seg, ori = cur
        if seg in visited:
            break
        visited.add(seg)
        path.append(cur)

        neighbors = [n for n in adj[cur] if n[0][0] in component]
        next_end, next_cigar = None, None
        for cand, cig in neighbors:
            if cand != prev:
                next_end = cand
                next_cigar = cig
                break

        if next_end is None:
            break

        cigars.append(next_cigar)
        prev = cur
        cur  = next_end

    return path, cigars

def merge_sequence(path, cigars):
    """
    Merge the path segments using simple <num>M CIGAR overlaps.
    """
    out = ""
    for i, (seg, ori) in enumerate(path):
        seq = segments[seg]
        if ori == "-":
            seq = revcomp(seq)
        if i == 0:
            out = seq
        else:
            m = re.match(r"(\d+)M", cigars[i-1])
            if not m:
                raise ValueError(f"Unsupported CIGAR: {cigars[i-1]}")
            overlap = int(m.group(1))
            out += seq[overlap:]
    return out

# ------------------------------------
# Find all connected components
# ------------------------------------

all_segments = set(segments.keys())
seen = set()
unitigs = []

for seg in all_segments:
    if seg in seen:
        continue

    comp = set()
    queue = deque([seg])
    while queue:
        s = queue.popleft()
        if s in comp:
            continue
        comp.add(s)
        for nbr in undirected[s]:
            if nbr not in comp:
                queue.append(nbr)

    seen |= comp

    # For this component, find the correct start end
    start_end = find_component_start(comp)
    path, cigars = walk_path(start_end, comp)

    if not path:
        continue

    seq = merge_sequence(path, cigars)
    unitigs.append(seq)

# ------------------------------------
# Emit FASTA
# ------------------------------------

for i, seq in enumerate(unitigs, 1):
    print(f">unitig_{i}")
    print(seq)
