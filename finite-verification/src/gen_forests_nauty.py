#!/usr/bin/env python3
"""Generate the subcubic carrier trees (|C| <= 9, Delta <= 3) for the finite
verification of Lemma (acyclic case) using nauty's gentreeg, the standard
trusted generator, instead of an in-house enumeration.

Requires the `gentreeg` binary (from the nauty suite) on PATH or in ./ .
Each tree is emitted by `gentreeg n -D3 -s` in sparse6 and parsed with
networkx. Substituting these into search_general.TREES and re-running the
search reproduces the zero-survivor result byte-for-byte (see paper, S:
Computational methodology).
"""
import subprocess, networkx as nx

def nauty_subcubic_trees(n, gentreeg="./gentreeg"):
    if n == 1:
        g = nx.Graph(); g.add_node(0); return [g]
    out = subprocess.run([gentreeg, str(n), "-D3", "-s", "-q"],
                         capture_output=True).stdout
    trees = []
    for line in out.split(b"\n"):
        line = line.strip()
        if not line or line.startswith(b">"):
            continue
        trees.append(nx.convert_node_labels_to_integers(nx.from_sparse6_bytes(line)))
    return trees

if __name__ == "__main__":
    for n in range(1, 10):
        print("n=%d : %d subcubic trees" % (n, len(nauty_subcubic_trees(n))))
