#!/usr/bin/python3


import numpy as np

from graphviz import Digraph
import networkx as nx
import networkx.algorithms as nxa
import re

import DeBruijnDNA
import AcyclicGraphDNA

import matplotlib.pyplot as plt

from collections import defaultdict




import argparse



parser = argparse.ArgumentParser()
parser.add_argument("input_file", type=str, help="Input file for the GAP problem")
parser.add_argument("output_file", type=str, help="Output file for the results")
parser.add_argument("solution_method", type=str, help="solve using")
parser.add_argument("load_file", type=str, help="load pkl file", default="{}")
parser.add_argument("sols_size", type=int, help="number of sols to evaluate", default=50)
parser.add_argument("fastas_file_number", type=int, help="fasta file number")
args = parser.parse_args()



def evaluate_penalty(qubo_dict, problem_size, *args):
    score = 0
    for i in range(problem_size):
        for j in range(i, problem_size):
            if i==j and args[0][i] == 1:
                score+=qubo_dict[(i,)]
            elif j>i and args[0][i] == 1 and args[0][j] == 1:
                score+=qubo_dict[(i, j)]
    return score


def compute_degrees(edge_list):
    in_deg = defaultdict(int)
    out_deg = defaultdict(int)
    nodes = set()

    for u, v in edge_list:
        out_deg[u] += 1
        in_deg[v] += 1
        nodes.add(u)
        nodes.add(v)

    # ensure zero-degree nodes appear
    for n in nodes:
        in_deg[n] += 0
        out_deg[n] += 0

    return in_deg, out_deg, nodes

def build_adj_out(edges):
    adj = defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
    return adj


def build_adj_in(edges):
    adj = defaultdict(list)
    for u, v in edges:
        adj[v].append(u)
    return adj


def longest_path_from(node, adj, node_lengths, memo, visiting):
    # Already solved
    res=0
    # print("current node", node_labels[node])
    if node in memo:
        return memo[node]

    # Cycle detected → cut
    if node in visiting:
        w = 0
        visiting = set()
        return res

    visiting.add(node)

    # Handle list-based node_lengths safely
    w = node_lengths[node]


    children = adj.get(node, [])


    if not children:
        res = res + w
    else:
        res = res + w + max(
            longest_path_from(v, adj, node_lengths, memo, visiting)
            for v in children
        )

    visiting = set()

    memo[node] = res
    return res

def enforce_linear(scored_edges):
    used_sinks = set()
    used_sources = set()
    kept = []

    for u, v in scored_edges:
        if u not in used_sinks and v not in used_sources:
            kept.append([u, v])
            used_sinks.add(u)
            used_sources.add(v)

    return kept

def write_solution_gfa(gfa_in, solve_edges, gfa_out):
    """
    gfa_in:   path to original GFA1 file (e.g., from miniasm)
    solve_edges: iterable of tuples like ('A+','B-'), oriented segment IDs
    gfa_out:  path to write the filtered GFA
    """
    s_lines = {}             # seg -> full S line
    l_lines = {}             # (from,fo,to,to_o) -> full L line

    # Parse just S and L
    with open(gfa_in) as fh:
        for line in fh:
            if line.startswith('S\t'):
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 3:
                    s_lines[parts[1]] = line.rstrip('\n')
            elif line.startswith('L\t'):
                p = line.rstrip('\n').split('\t')
                if len(p) >= 6:
                    l_lines[(p[1], p[2], p[3], p[4])] = line.rstrip('\n')

    # Keep only edges in solve_edges (must match orientation exactly)
    keep_segments = set()
    keep_links = []
    for u, v in solve_edges:
        us, uo = u[:-1], u[-1]   # e.g. 'A+', -> ('A','+')
        vs, vo = v[:-1], v[-1]
        key = (us, uo, vs, vo)
        if key in l_lines:
            keep_links.append(l_lines[key])
            keep_segments.update([us, vs])

    # Write minimal GFA: header, needed S-lines, then chosen L-lines
    with open(gfa_out, 'w') as out:
        out.write('H\tVN:Z:1.0\n')
        for seg in sorted(keep_segments):
            out.write((s_lines.get(seg, f"S\t{seg}\t*\n")).rstrip('\n') + '\n')
        for l in keep_links:
            out.write(l.rstrip('\n') + '\n')

def longest_nodes_from(node, adj, memo, visiting):
    # Already solved
    res=0
    # print("current node", node_labels[node])
    if node in memo:
        return memo[node]

    # Cycle detected → cut
    if node in visiting:
        w = 0
        visiting = set()
        return res

    visiting.add(node)

    # Handle list-based node_lengths safely
    w = 1


    children = adj.get(node, [])


    if not children:
        res = res + w
    else:
        res = res + w + max(
            longest_nodes_from(v, adj, memo, visiting)
            for v in children
        )

    visiting = set()

    memo[node] = res
    return res



filename = args.input_file
# filename = "graph.gfa"
segments, links, containments = AcyclicGraphDNA.load_file(filename)
segments = [l.replace(":", "/", 1) for l in segments]
links = [l.replace(":", "/", 2) for l in links]


initial_graph, adjacency_matrix, nodes_indices, nodes_labels, links_edges = AcyclicGraphDNA.get_initial_graph(segments, links)
strand_graph = AcyclicGraphDNA.get_strand_graph(initial_graph)

strand_graph_adjacency_matrix = nx.to_numpy_matrix(strand_graph)

edges = np.argwhere(strand_graph_adjacency_matrix==1).tolist()
nodes = list(strand_graph.nodes)

strand_graph_undirected = AcyclicGraphDNA.complement_to_undirected_graph(strand_graph)
adjacency_matrix = nx.to_numpy_matrix(strand_graph_undirected)
part_map, _ , _ = AcyclicGraphDNA.part_graph(cluster_number=1, adjacency_matrix=adjacency_matrix)

graphs = []
node_labels = list(strand_graph_undirected.nodes)
node_lengths = [len(segments[i//2])-len(node_labels[i])-3 for i in range(len(node_labels))]
for part, nodes_indices in part_map.items():
    labels = [node_labels[i] for i in nodes_indices]
    graphs.append(strand_graph.subgraph(labels))



if args.solution_method == "load":
    import dill as pickle
    with open(args.load_file, "rb") as f:
        sols = pickle.load(f)
    for pg in graphs:
        adjacency_matrix = nx.to_numpy_matrix(pg)
        edges_indices = np.argwhere(adjacency_matrix==1).tolist()

        Q = AcyclicGraphDNA.to_qubo(adjacency_matrix, A=2)
    qubo_dict = defaultdict(int)
    for i in range(2313):
        for j in range(i, 2313):
            if i==j:
                qubo_dict[(i,)] = Q[i][j]
            else:
                qubo_dict[(i, j)] = Q[i][j]
    avg_sol = [evaluate_penalty(qubo_dict, 2313, sols[2][i]) for i in range(args.sols_size)]
    average = np.mean(avg_sol)
    best = np.min(avg_sol)
    probsols = evaluate_penalty(qubo_dict, 2313, sols[0])

i = args.fastas_file_number
print("Processing solution ", i)

if args.solution_method == "dwave":
    sol_no = 2*i
    solve_edges=[]
    spins, energy = [solution[sol_no][i] for i in solution[sol_no].keys()], solution[sol_no+1]

    solve_edges_idices=[]
    for key in range(len(spins)):
        if spins[key]==1:
            solve_edges_idices.append(edges_indices[key])

    node_labels = list(pg.nodes)
    solve_edges_part = [(node_labels[e[0]], node_labels[e[1]]) for e in solve_edges_idices]
    solve_edges+=solve_edges_part
elif args.solution_method == "HADOF" or args.solution_method == "load":
    if i == 0:
        solve_edges=[]
        solve_edges_idices=[]
        for key in range(len(sols[0])):
            if sols[0][key]==1:
                solve_edges_idices.append(edges_indices[key])

        node_labels = list(pg.nodes)
        solve_edges_part = [(node_labels[e[0]], node_labels[e[1]]) for e in solve_edges_idices]
        solve_edges+=solve_edges_part
    else:
        sol_no = i-1
        solve_edges=[]
        solve_edges_idices=[]
        for key in range(len(sols[2][sol_no])):
            if sols[2][sol_no][key]==1:
                solve_edges_idices.append(edges_indices[key])

        node_labels = list(pg.nodes)
        solve_edges_part = [(node_labels[e[0]], node_labels[e[1]]) for e in solve_edges_idices]
        solve_edges+=solve_edges_part



in_deg, out_deg, sub_nodes = compute_degrees(solve_edges_idices)
check_nodes_out = [n for n in sub_nodes if out_deg[n] > 1]
check_edges_out = [[u, v] for u, v in solve_edges_idices if u in check_nodes_out]
check_nodes_in = [n for n in sub_nodes if in_deg[n] > 1]
check_edges_in = [[u, v] for u, v in solve_edges_idices if v in check_nodes_in]

adj_out = build_adj_out(solve_edges_idices)
adj_in = build_adj_in(solve_edges_idices)


memo = {}

targets = check_edges_out

path_len = {
    v: longest_path_from(v, adj_out, node_lengths, memo, set())
    for _, v in targets
}

scored = sorted(
    check_edges_out,
    key=lambda e: path_len[e[1]],
    reverse=True
)

recovered_edges_out = enforce_linear(scored)

memo = {}

targets = check_edges_in

path_len = {
    u: longest_path_from(u, adj_in, node_lengths, memo, set())
    for u, _ in targets
}

scored = sorted(
    check_edges_in,
    key=lambda e: path_len[e[0]],
    reverse=True
)
recovered_edges_in = enforce_linear(scored)
repaired_edges_out = [[u, v] for u, v in solve_edges_idices if ([u, v] not in check_edges_out or [u, v] in recovered_edges_out)]
repaired_edges = [[u, v] for u, v in repaired_edges_out if [u, v] not in check_edges_in or [u, v] in recovered_edges_in]

in_deg, out_deg, sub_nodes = compute_degrees(repaired_edges)
sinks   = {n for n in sub_nodes if out_deg[n] == 0}
sources = {n for n in sub_nodes if in_deg[n] == 0}
bad_nodes = {n for n in nodes_indices if n not in sub_nodes}
sinks = sinks.union(bad_nodes)
sources = sources.union(bad_nodes)

reconnecting_edges = [
    [u, v]
    for u, v in edges
    if u in sinks and v in sources
]

adj_out = build_adj_out(repaired_edges)
# adj_out = build_adj_out(solve_edges_idices)
memo = {}
path_len = {v: longest_path_from(v, adj_out, node_lengths, memo, set()) for _, v in reconnecting_edges}
scored = sorted(
    reconnecting_edges,
    key=lambda e: path_len[e[1]],
    reverse=True
)

connections = enforce_linear(scored)

final_edges = repaired_edges + connections

in_deg, out_deg, sub_nodes = compute_degrees(final_edges)
sinks   = {n for n in sub_nodes if out_deg[n] == 0}
sources = {n for n in sub_nodes if in_deg[n] == 0}

memo={}
adj_out = build_adj_out(final_edges)
path_len = {v: longest_nodes_from(v, adj_out, memo, set()) for v in sources}

# chosen_paths = [list(path_len.keys())[i] for i in range(len(path_len)) if list(path_len.values())[i]>40] 
chosen_paths = [max(path_len, key=path_len.get)]
if chosen_paths:

    newfinal_edges = []
    for path in chosen_paths:
        node = path
        while node in adj_out and adj_out[node]:
            newfinal_edges.append([node, adj_out[node][0]])
            node = adj_out[node][0]   

    solve_edges = []
    node_labels = list(pg.nodes)
    solve_edges_part = [(node_labels[e[0]], node_labels[e[1]]) for e in newfinal_edges]
    solve_edges+=solve_edges_part

    out_edges = [(u.replace('/', ':'), v.replace('/', ':')) for (u, v) in solve_edges]
    write_solution_gfa(filename, out_edges, args.output_file+"_"+str(i)+".gfa")


