#!/usr/bin/env python3
"""Generalized search: forest-carrier case for deficiency d in {1,2,3} and
interface budget c <= 5 - d.  R-attachments enumerated definition-level:
every rho in R with C-edges has an alpha-invariant C-neighbourhood (any
multiset of invariant sets, sizes 1..3, total <= cmax)."""
import sys
from itertools import combinations, combinations_with_replacement, product
from collections import Counter
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher

def partitions(n, maxpart=None):
    maxpart = maxpart if maxpart is not None else n
    if n == 0:
        yield []; return
    for p in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - p, p):
            yield [p] + rest

def trees_of_size(s):
    if s == 1:
        g = nx.Graph(); g.add_node(0); return [g]
    if s == 2:
        g = nx.Graph(); g.add_edge(0, 1); return [g]
    return list(nx.nonisomorphic_trees(s))

TREES = {s: [t for t in trees_of_size(s)
             if max((d for _, d in t.degree()), default=0) <= 3]
         for s in range(1, 11)}

def forests(m, budget):
    for part in partitions(m):
        k = len(part)
        if m + 2 * k > budget: continue
        cnt = sorted(Counter(part).items())
        cls = [(s, list(combinations_with_replacement(range(len(TREES[s])), t)))
               for s, t in cnt]
        for combo in product(*[o for _, o in cls]):
            G = nx.Graph(); off = 0
            for (s, _), idxs in zip(cls, combo):
                for i in idxs:
                    t = TREES[s][i]
                    G.add_nodes_from(range(off, off + s))
                    G.add_edges_from((u + off, v + off) for u, v in t.edges())
                    off += s
            yield G, k

def nontrivial_autos(T):
    for mp in GraphMatcher(T, T).isomorphisms_iter():
        if any(mp[v] != v for v in T.nodes()):
            yield dict(mp)

def invariant_sets(T, alpha):
    out = []
    nodes = list(T.nodes())
    for sz in (1, 2, 3):
        for S in combinations(nodes, sz):
            fs_ = frozenset(S)
            if frozenset(alpha[x] for x in S) == fs_:
                out.append(fs_)
    return out

def r_collections(invsets, cmax):
    res = [[]]
    def rec(start, cur, tot):
        for i in range(start, len(invsets)):
            S = invsets[i]
            if tot + len(S) <= cmax:
                nxt = cur + [S]
                res.append(nxt)
                rec(i, nxt, tot + len(S))
    rec(0, [], 0)
    return res

def girth_ok(adj):
    nodes = list(adj)
    for u, v in combinations(nodes, 2):
        common = adj[u] & adj[v]
        if v in adj[u] and common: return False
        if len(common) >= 2: return False
    return True

def cycles56(adj):
    nodes = sorted(adj, key=str)
    idx = {v: i for i, v in enumerate(nodes)}
    out = []
    for s in nodes:
        stack = [(s, (s,))]
        while stack:
            v, path = stack.pop()
            for w in adj[v]:
                if w == s and len(path) >= 3:
                    if len(path) in (3, 4):
                        raise AssertionError('short cycle')
                    if idx[path[1]] < idx[path[-1]]:
                        out.append(path)
                elif idx.get(w, -1) > idx[s] and w not in path and len(path) < 6:
                    stack.append((w, path + (w,)))
    return out

def ipr_verdict(adj, C_nodes):
    def complete(v):
        return (v in C_nodes) or (str(v).startswith('d') and len(adj[v]) == 3)
    cyc = cycles56(adj)
    pents = [c for c in cyc if len(c) == 5]
    for a, b in combinations(pents, 2):
        if set(a) & set(b): return None
    cert_hex = []
    for c in cyc:
        if len(c) != 6: continue
        ok = True
        for i, j in combinations(range(6), 2):
            if j - i in (1, 5): continue
            u, v = c[i], c[j]
            if v in adj[u]: raise AssertionError('chord')
            if not (complete(u) or complete(v)):
                ok = False; break
        if ok: cert_hex.append(c)
    faces = [('P', c) for c in pents] + [('H', c) for c in cert_hex]
    def es(c): return {frozenset((c[i], c[(i+1) % len(c)])) for i in range(len(c))}
    for (ta, a), (tb, b) in combinations(faces, 2):
        I = set(a) & set(b)
        if not I: continue
        E = es(a) & es(b)
        if not (len(I) == 2 and len(E) == 1 and set(next(iter(E))) == I): return None
        if ta == 'P' and tb == 'P': return None
    vc, ec = Counter(), Counter()
    for _, c in faces:
        for v in c: vc[v] += 1
        for e in es(c): ec[e] += 1
    if vc and max(vc.values()) > 3: return None
    if ec and max(ec.values()) > 2: return None
    return True

def run(d, cmax):
    DELTAS = ['d%d' % i for i in range(d)]
    stats = Counter(); survivors = []
    budget = cmax + 3 * d
    for m in range(2, budget - 2 + 1):
        for T, k in forests(m, budget):
            degT = dict(T.degree())
            for alpha in nontrivial_autos(T):
                stats['(T,alpha)'] += 1
                M = {v for v in T if alpha[v] != v}
                FixC = set(T) - M
                invs = invariant_sets(T, alpha)
                for ratt in r_collections(invs, cmax):
                    stats['(T,alpha,r)'] += 1
                    rdeg = Counter()
                    for S in ratt:
                        for x in S: rdeg[x] += 1
                    if any(degT[x] + rdeg[x] > 3 for x in T):
                        stats['kill deg'] += 1; continue
                    bad = False
                    for x in M:
                        tot = len(set(T[x]) & FixC) + rdeg[x]
                        if tot > 1 or (alpha[x] in T[x] and tot > 0):
                            bad = True; break
                    if bad:
                        stats['kill anchor'] += 1; continue
                    s = {x: 3 - degT[x] - rdeg[x] for x in T}
                    if sum(s.values()) > 3 * d:
                        stats['kill S_D budget'] += 1; continue
                    adjB = {v: set(T[v]) for v in T}
                    for i, S in enumerate(ratt):
                        lbl = 'r%d' % i
                        adjB[lbl] = set(S)
                        for x in S: adjB[x].add(lbl)
                    if not girth_ok(adjB):
                        stats['kill girth base'] += 1; continue
                    Bx = nx.Graph()
                    for v, nb in adjB.items():
                        Bx.add_node(v); Bx.add_edges_from((v, w) for w in nb)
                    dist = dict(nx.all_pairs_shortest_path_length(Bx))
                    verts = sorted(T, key=lambda v: -s[v])
                    dsets = {dl: set() for dl in DELTAS}
                    pair_used = set()
                    def bt(i):
                        if i == len(verts):
                            finish(); return
                        x = verts[i]
                        for chosen in combinations(DELTAS, s[x]):
                            ok = True; newpairs = []
                            for dl in chosen:
                                if len(dsets[dl]) >= 3: ok = False; break
                                for y in dsets[dl]:
                                    if dist[x].get(y, 99) < 3: ok = False; break
                                    fp = frozenset((x, y))
                                    if fp in pair_used or fp in newpairs: ok = False; break
                                    newpairs.append(fp)
                                if not ok: break
                            if not ok: continue
                            for dl in chosen: dsets[dl].add(x)
                            for fp in newpairs: pair_used.add(fp)
                            bt(i + 1)
                            for dl in chosen: dsets[dl].remove(x)
                            for fp in newpairs: pair_used.discard(fp)
                    def finish():
                        stats['complete (pre-IPR)'] += 1
                        adjK = {v: set(nb) for v, nb in adjB.items()}
                        for dl in DELTAS:
                            if dsets[dl]:
                                adjK[dl] = set(dsets[dl])
                                for x in dsets[dl]: adjK[x].add(dl)
                        if ipr_verdict(adjK, set(T.nodes())) is None:
                            stats['kill IPR faces'] += 1
                        else:
                            stats['SURVIVORS'] += 1
                            survivors.append((m, sorted(T.edges()), dict(alpha), ratt,
                                              {dl: sorted(dsets[dl]) for dl in DELTAS}))
                    bt(0)
    print('--- d=%d, c<=%d (|C| <= %d) ---' % (d, cmax, budget - 2))
    for key in sorted(stats): print('  %-28s %d' % (key, stats[key]))
    if survivors:
        print('  !!! survivors:')
        for sv in survivors: print('   ', sv)
    else:
        print('  no survivors: impossible.')
    return survivors

if __name__ == '__main__':
    total = 0
    for d, cmax in [(3, 2), (2, 3), (1, 4)]:
        total += len(run(d, cmax))
    print()
    print('TOTAL SURVIVORS over all regimes:', total)
