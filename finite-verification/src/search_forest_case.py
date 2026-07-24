#!/usr/bin/env python3
"""
Exhaustive configuration search for the forest-carrier case of the unified
hidden-symmetry theorem for asymmetric IPR fullerenes.

Setting.  F is an IPR fullerene: cubic, planar, 3-connected, girth 5, and
(Kardos-Skrekovski) every cyclic 5- or 6-edge cut is trivial.  Consequences
used as filters (each cited/proved in the accompanying write-up):
  (S1) girth 5: no 3- or 4-cycles;
  (S2) every 5-cycle bounds a (pentagonal) face;
  (S3) two faces sharing a vertex share an edge (cubic); two pentagon faces
       sharing an edge violate IPR; hence two distinct 5-cycles sharing even
       one vertex are impossible;
  (S4) every chordless 6-cycle bounds a (hexagonal) face;
  (S5) two distinct faces intersect in nothing or exactly one edge (with its
       two endpoints); each vertex lies on exactly 3 faces, each edge on 2.

Configuration.  phi = alpha  cup  id_R with dom(phi) = C ⊔ R, alpha a
nontrivial automorphism of F[C], R fixed pointwise, F[C] a forest,
c = e(C,R) <= 2, S_D = V \\ dom(phi), |S_D| = d(F) = d in {2,3}.
The run with d = 3 (three deleted vertices delta0,delta1,delta2, each of
degree 3) subsumes d = 2, since any d = 2 configuration is a d = 3
configuration that happens not to use the third deleted vertex.

Known counting facts encoded:
  * every x in C has degree exactly 3 in F: deg_T(x) + r(x) + s(x) = 3;
  * |C| + 2k = c + e(C,S_D) <= 2 + 9 = 11  (k = #components of the forest),
    hence |C| <= 9;
  * a fixed neighbour of a moved vertex x is a common neighbour of x and
    alpha(x) (anchor lemma); with girth 5 a moved vertex has at most ONE
    fixed neighbour, and none if x ~ alpha(x);
  * an R-vertex rho is fixed, so N_C(rho) is alpha-invariant; in particular
    an R-edge to a moved vertex forces rho to be adjacent to the whole
    orbit, which must therefore be a 2-orbit and consume the entire budget
    c = 2.

The search enumerates all (forest T, nontrivial alpha in Aut(T), R-attachment,
S_D-attachment with identifications) and reports every configuration that
survives all sound filters.  Zero survivors = the forest case is impossible.
"""

import sys
from itertools import combinations, combinations_with_replacement, product
from collections import Counter, defaultdict
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher

D = 3                # |S_D| (d = 3 run subsumes d = 2)
DELTAS = ['d0', 'd1', 'd2']
DELTA_CAP = 3        # each deleted vertex has degree 3
C_BUDGET = 2         # c = e(C,R) <= 2

# ---------------------------------------------------------------- forests

def partitions(n, maxpart=None):
    maxpart = maxpart if maxpart is not None else n
    if n == 0:
        yield []
        return
    for p in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - p, p):
            yield [p] + rest

def trees_of_size(s):
    if s == 1:
        g = nx.Graph(); g.add_node(0); return [g]
    if s == 2:
        g = nx.Graph(); g.add_edge(0, 1); return [g]
    return list(nx.nonisomorphic_trees(s))

TREES = {}
for s in range(1, 10):
    TREES[s] = [t for t in trees_of_size(s)
                if max((d for _, d in t.degree()), default=0) <= 3]

def forests(m):
    """All forests on m vertices, max degree <= 3, up to isomorphism,
    satisfying m + 2k <= 11 (else e(C,S_D) budget already violated)."""
    for part in partitions(m):
        k = len(part)
        if m + 2 * k > C_BUDGET + 3 * D:      # m + 2k <= c + 9 <= 11
            continue
        cnt = sorted(Counter(part).items())
        choice_lists = []
        for size, mult in cnt:
            opts = list(combinations_with_replacement(range(len(TREES[size])), mult))
            choice_lists.append((size, opts))
        for combo in product(*[opts for _, opts in choice_lists]):
            G = nx.Graph()
            off = 0
            for (size, _), idxs in zip(choice_lists, combo):
                for i in idxs:
                    t = TREES[size][i]
                    G.add_nodes_from(range(off, off + size))
                    G.add_edges_from((u + off, v + off) for u, v in t.edges())
                    off += size
            yield G, k

# ------------------------------------------------------- automorphisms

def nontrivial_autos(T):
    gm = GraphMatcher(T, T)
    for mp in gm.isomorphisms_iter():
        if any(mp[v] != v for v in T.nodes()):
            yield dict(mp)

# ------------------------------------------------------- R attachments

def r_options(T, alpha, M, FixC):
    """Yield lists of (rho_label, tuple_of_C_endpoints); total edges <= 2.
    Single R-edges may only go to fixed vertices (an R-edge to a moved x
    forces rho ~ whole orbit of x: needs orbit length 2 and both edges)."""
    opts = [[]]
    fixcap = [x for x in FixC if T.degree(x) <= 2]
    # c = 1
    for x in fixcap:
        opts.append([('r0', (x,))])
    # c = 2, two distinct rho's, endpoints fixed
    for x, y in combinations_with_replacement(fixcap, 2):
        if x == y:
            if T.degree(x) <= 1:
                opts.append([('r0', (x,)), ('r1', (x,))])
        else:
            opts.append([('r0', (x,)), ('r1', (y,))])
    # c = 2, one rho with two fixed endpoints
    for x, y in combinations(fixcap, 2):
        opts.append([('r0', (x, y))])
    # c = 2, one rho covering a 2-orbit of alpha
    seen = set()
    for x in M:
        y = alpha[x]
        if y != x and alpha[y] == x and not T.has_edge(x, y) \
           and T.degree(x) <= 2 and T.degree(y) <= 2:
            key = frozenset((x, y))
            if key not in seen:
                seen.add(key)
                opts.append([('r0', (x, y))])
    return opts

# ------------------------------------------------------------ girth check

def girth_ok(adj):
    """No 3- or 4-cycles among known edges: no adjacent pair with a common
    neighbour, no pair with two common neighbours."""
    nodes = list(adj)
    for u, v in combinations(nodes, 2):
        common = adj[u] & adj[v]
        if v in adj[u] and common:
            return False            # triangle
        if len(common) >= 2:
            return False            # 4-cycle
    return True

# ------------------------------------------------------- cycle enumeration

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
                        raise AssertionError('short cycle escaped girth filter')
                    if idx[path[1]] < idx[path[-1]]:
                        out.append(path)
                elif idx.get(w, -1) > idx[s] and w not in path and len(path) < 6:
                    stack.append((w, path + (w,)))
    return out

# --------------------------------------------------------- IPR face logic

def ipr_verdict(adj, C_nodes):
    """Return None if configuration is contradictory, else a dict of face data
    (a genuine survivor)."""
    def complete(v):
        return (v in C_nodes) or (str(v).startswith('d') and len(adj[v]) == 3)

    cyc = cycles56(adj)
    pents = [c for c in cyc if len(c) == 5]
    hexes6 = [c for c in cyc if len(c) == 6]

    # (S2)+(S3): two distinct 5-cycles sharing a vertex -> adjacent pentagons
    for a, b in combinations(pents, 2):
        if set(a) & set(b):
            return None

    # (S4): certified chordless 6-cycles are hexagonal faces
    cert_hex = []
    for c in hexes6:
        chordless_certified = True
        n = len(c)
        for i, j in combinations(range(n), 2):
            if j - i in (1, n - 1):
                continue                       # cycle edge
            u, v = c[i], c[j]
            if v in adj[u]:
                raise AssertionError('chord present: girth filter should have fired')
            if not (complete(u) or complete(v)):
                chordless_certified = False
                break
        if chordless_certified:
            cert_hex.append(c)

    faces = [('P', c) for c in pents] + [('H', c) for c in cert_hex]

    def edgeset(c):
        return {frozenset((c[i], c[(i + 1) % len(c)])) for i in range(len(c))}

    # (S5) + (S3): pairwise face intersections
    for (ta, a), (tb, b) in combinations(faces, 2):
        I = set(a) & set(b)
        if not I:
            continue
        E = edgeset(a) & edgeset(b)
        ok = (len(I) == 2 and len(E) == 1 and set(next(iter(E))) == I)
        if not ok:
            return None                        # faces meet illegally
        if ta == 'P' and tb == 'P':
            return None                        # adjacent pentagons (IPR)

    # vertex on <= 3 faces, edge on <= 2 faces
    vcount, ecount = Counter(), Counter()
    for _, c in faces:
        for v in c:
            vcount[v] += 1
        for e in edgeset(c):
            ecount[e] += 1
    if vcount and max(vcount.values()) > 3:
        return None
    if ecount and max(ecount.values()) > 2:
        return None

    return {'pentagons': pents, 'hexagons': cert_hex}

# ----------------------------------------------------------------- search

def run():
    stats = Counter()
    survivors = []

    for m in range(2, 10):
        for T, k in forests(m):
            degT = dict(T.degree())
            for alpha in nontrivial_autos(T):
                stats['(T,alpha) pairs'] += 1
                M = {v for v in T if alpha[v] != v}
                FixC = set(T) - M
                # anchor pre-check on T alone
                bad = False
                for x in M:
                    ft = len(set(T[x]) & FixC)
                    if ft > 1 or (alpha[x] in T[x] and ft > 0):
                        bad = True
                        break
                if bad:
                    stats['killed: anchor (T only)'] += 1
                    continue

                for ropt in r_options(T, alpha, M, FixC):
                    stats['(T,alpha,r) triples'] += 1
                    rdeg = Counter()
                    for _, ends in ropt:
                        for x in ends:
                            rdeg[x] += 1
                    # anchor rule including R-edges
                    bad = False
                    for x in M:
                        tot = len(set(T[x]) & FixC) + rdeg[x]
                        if tot > 1 or (alpha[x] in T[x] and tot > 0):
                            bad = True
                            break
                    if bad:
                        stats['killed: anchor (with R)'] += 1
                        continue
                    s = {x: 3 - degT[x] - rdeg[x] for x in T}
                    if any(v < 0 for v in s.values()):
                        stats['killed: degree overflow'] += 1
                        continue
                    if sum(s.values()) > D * DELTA_CAP:
                        stats['killed: S_D budget'] += 1
                        continue

                    # base graph B = T + R-edges
                    adjB = {v: set(T[v]) for v in T}
                    for rho, ends in ropt:
                        adjB[rho] = set(ends)
                        for x in ends:
                            adjB[x].add(rho)
                    if not girth_ok(adjB):
                        stats['killed: girth (base)'] += 1
                        continue

                    # distances in base graph
                    Bx = nx.Graph()
                    for v, nb in adjB.items():
                        Bx.add_node(v)
                        Bx.add_edges_from((v, w) for w in nb)
                    dist = dict(nx.all_pairs_shortest_path_length(Bx))

                    # backtracking over delta assignments
                    verts = sorted(T, key=lambda v: -s[v])
                    dsets = {dl: set() for dl in DELTAS}
                    pair_used = set()

                    def bt(i):
                        if i == len(verts):
                            finish()
                            return
                        x = verts[i]
                        need = s[x]
                        # no symmetry pruning: enumerate all label subsets
                        # (label-permutation duplicates are harmless)
                        for chosen in combinations(DELTAS, need):
                            ok = True
                            newpairs = []
                            for dl in chosen:
                                if len(dsets[dl]) >= DELTA_CAP:
                                    ok = False
                                    break
                                for y in dsets[dl]:
                                    if dist[x].get(y, 99) < 3:
                                        ok = False
                                        break
                                    fp = frozenset((x, y))
                                    if fp in pair_used or fp in newpairs:
                                        ok = False
                                        break
                                    newpairs.append(fp)
                                if not ok:
                                    break
                            if not ok:
                                continue
                            for dl in chosen:
                                dsets[dl].add(x)
                            for fp in newpairs:
                                pair_used.add(fp)
                            bt(i + 1)
                            for dl in chosen:
                                dsets[dl].remove(x)
                            for fp in newpairs:
                                pair_used.discard(fp)

                    def finish():
                        stats['complete assignments (pre-IPR survivors)'] += 1
                        adjK = {v: set(nb) for v, nb in adjB.items()}
                        for dl in DELTAS:
                            if dsets[dl]:
                                adjK[dl] = set(dsets[dl])
                                for x in dsets[dl]:
                                    adjK[x].add(dl)
                        assert girth_ok(adjK)
                        verdict = ipr_verdict(adjK, set(T.nodes()))
                        if verdict is None:
                            stats['killed: IPR/planarity faces'] += 1
                        else:
                            stats['SURVIVORS'] += 1
                            survivors.append({
                                'm': m, 'k': k,
                                'T_edges': sorted(map(tuple, T.edges())),
                                'alpha': {v: alpha[v] for v in T if alpha[v] != v},
                                'r': ropt,
                                'deltas': {dl: sorted(dsets[dl]) for dl in DELTAS if dsets[dl]},
                                'faces': verdict,
                            })

                    bt(0)

    print('=== exhaustive forest-carrier search, d = %d (subsumes d = 2), c <= %d ===' % (D, C_BUDGET))
    for key in sorted(stats):
        print('%-45s %d' % (key, stats[key]))
    print()
    if survivors:
        print('!!! %d surviving configuration(s):' % len(survivors))
        for sv in survivors:
            print(sv)
    else:
        print('No configuration survives: the forest-carrier case is impossible. QED')
    return survivors

if __name__ == '__main__':
    run()
