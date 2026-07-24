#!/usr/bin/env python3
"""Construct a cubic, 3-connected, girth-5 graph carrying a localised
nontrivial partial automorphism of deficiency 3 behind a 2-edge interface.
Fragment: path x-u-w-v-y, psi = reversal (x<->y, u<->v, w fixed);
rho ~ x,y (fixed, the 2-edge interface); deleted s1~u,y  s2~v,x  s3~w."""
import random, networkx as nx
from itertools import combinations

C = ['x','u','w','v','y']
FRAG = [('x','u'),('u','w'),('w','v'),('v','y'),
        ('rho','x'),('rho','y'),
        ('s1','u'),('s1','y'),('s2','v'),('s2','x'),('s3','w')]
PSI = {'x':'y','y':'x','u':'v','v':'u','w':'w'}
DEL = {'s1','s2','s3'}

def is_pa(G, phi):
    dom = set(phi)
    for a, b in combinations(sorted(dom), 2):
        if G.has_edge(a, b) != G.has_edge(phi[a], phi[b]):
            return False
    return True

def try_build(m, rng, tries=4000):
    extra = ['e%d' % i for i in range(m)]
    need0 = {'rho':1, 's1':1, 's2':1, 's3':2}
    for _ in range(tries):
        G = nx.Graph(); G.add_edges_from(FRAG); G.add_nodes_from(extra)
        stubs = []
        for v, k in need0.items(): stubs += [v] * k
        stubs += [v for v in extra for _ in range(3)]
        rng.shuffle(stubs)
        ok = True
        # greedy matching respecting simplicity + girth 5
        while stubs:
            a = stubs.pop()
            cand = [i for i, b in enumerate(stubs)
                    if b != a and not G.has_edge(a, b)
                    and (nx.shortest_path_length(G, a, b) >= 4
                         if nx.has_path(G, a, b) else True)]
            if not cand: ok = False; break
            i = rng.choice(cand); b = stubs.pop(i)
            G.add_edge(a, b)
        if not ok: continue
        if any(d != 3 for _, d in G.degree()): continue
        if nx.girth(G) < 5: continue
        if not nx.is_connected(G) or nx.node_connectivity(G) < 3: continue
        phi = {v: PSI[v] for v in C}
        for v in G.nodes():
            if v not in C and v not in DEL: phi[v] = v
        if not is_pa(G, phi): continue
        return G, phi
    return None, None

rng = random.Random(20260721)
for m in range(3, 16):
    if (5 + 3 * m) % 2: continue   # parity: total stub count must be even
    G, phi = try_build(m, rng)
    if G:
        n = G.number_of_nodes()
        moved = [v for v in phi if phi[v] != v]
        print("FOUND witness: n =", n, "| girth =", nx.girth(G),
              "| connectivity =", nx.node_connectivity(G))
        print("  deficiency k =", n - len(phi), "| |supp| =", len(moved))
        print("  interface edges C-R:",
              sum(1 for a, b in G.edges() if (a in C) != (b in C)
                  and (a not in DEL and b not in DEL)))
        print("  |Aut(G)| =", sum(1 for _ in
              nx.algorithms.isomorphism.GraphMatcher(G, G).isomorphisms_iter()))
        print("  edges =", sorted(tuple(sorted(e)) for e in G.edges()))
        # confirm the two 5-cycles sharing two edges
        p1 = ['x','u','w','v','s2']; p2 = ['u','w','v','y','s1']
        def isc(c): return all(G.has_edge(c[i], c[(i+1) % len(c)]) for i in range(len(c)))
        print("  5-cycle x-u-w-v-s2:", isc(p1), "| 5-cycle u-w-v-y-s1:", isc(p2),
              "| shared edges:", 2)
        break
