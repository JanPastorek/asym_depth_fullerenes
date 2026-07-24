#!/usr/bin/env python3
"""Verify the explicit 20-vertex witness graph G of the paper's
Proposition (girth-5 witness), exactly as listed in the Computational
methodology section.

G is claimed to be: cubic, 3-connected, girth 5, with trivial automorphism
group, non-planar, and to admit the partial automorphism that fixes every
vertex outside {x,u,v,y} u {s1,s2,s3} and acts as psi = (x y)(u v); this
partial automorphism has rank 17, support {x,u,v,y}, and the only edges
between C = {x,u,w,v,y} and the pointwise-fixed remainder are rho-x and rho-y.
"""
from itertools import combinations
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher

EDGES = [
    ("x","u"),("u","w"),("w","v"),("v","y"),("rho","x"),("rho","y"),
    ("s1","u"),("s1","y"),("s2","v"),("s2","x"),("s3","w"),
    ("e0","e3"),("e0","e4"),("e0","s3"),("e1","e4"),("e1","e7"),("e1","s2"),
    ("e2","e3"),("e2","e7"),("e2","e9"),("e3","rho"),
    ("e4","e5"),("e5","e6"),("e5","e9"),("e6","e10"),("e6","s1"),
    ("e7","e10"),("e8","e9"),("e8","e10"),("e8","s3"),
]
DELETED = {"s1","s2","s3"}
PSI = {"x":"y","y":"x","u":"v","v":"u"}      # on moved vertices; w is fixed

def partial_aut(G):
    """phi: fix all of V\\DELETED, applying PSI where defined."""
    dom = set(G) - DELETED
    phi = {z: PSI.get(z, z) for z in dom}
    return phi

def is_partial_automorphism(G, phi):
    dom = sorted(phi)
    for a, b in combinations(dom, 2):
        if G.has_edge(a, b) != G.has_edge(phi[a], phi[b]):
            return False
    return True

def aut_count(G, cap=2):
    n = 0
    for _ in GraphMatcher(G, G).isomorphisms_iter():
        n += 1
        if n >= cap:
            break
    return n

def main():
    G = nx.Graph(); G.add_edges_from(EDGES)
    checks = {}
    checks["n = 20"]                 = (G.number_of_nodes() == 20)
    checks["cubic (all deg 3)"]      = all(d == 3 for _, d in G.degree())
    checks["3-connected"]            = (nx.node_connectivity(G) == 3)
    checks["girth 5"]                = (nx.girth(G) == 5)
    checks["trivial Aut (asymmetric)"] = (aut_count(G) == 1)
    checks["non-planar"]             = (not nx.check_planarity(G)[0])

    phi = partial_aut(G)
    checks["phi is a partial automorphism"] = is_partial_automorphism(G, phi)
    checks["rank(phi) = 17"]         = (len(phi) == 17)
    supp = {z for z in phi if phi[z] != z}
    checks["support = {x,u,v,y}"]    = (supp == {"x","u","v","y"})

    C = {"x","u","w","v","y"}
    fixed_remainder = (set(G) - DELETED) - C
    interface = [(a, b) for a, b in G.edges()
                 if (a in C) != (b in C)
                 and (a in fixed_remainder or b in fixed_remainder)]
    checks["interface C--fixed = {rho-x, rho-y}"] = (
        {frozenset(e) for e in interface} == {frozenset(("rho","x")),
                                              frozenset(("rho","y"))})

    width = max(len(k) for k in checks)
    ok = True
    for k, v in checks.items():
        print(f"  [{'OK' if v else 'FAIL'}] {k}")
        ok = ok and v
    print("\nALL CHECKS PASSED" if ok else "\nSOME CHECK FAILED")
    return 0 if ok else 1

if __name__ == "__main__":
    raise SystemExit(main())
