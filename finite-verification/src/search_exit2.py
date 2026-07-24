#!/usr/bin/env python3
"""Clause (2') search, v2 (fast slot backtracking).  See search_exit.py header.
Usage: search_exit2.py d cmax mlo mhi [mod rem]"""
import sys
from itertools import combinations, product
from collections import Counter, defaultdict
import networkx as nx
from search_general import forests, girth_ok, ipr_verdict

def collections_of(cands, cmax):
    res = [[]]
    def rec(start, cur, tot):
        for i in range(start, len(cands)):
            S = cands[i]
            if tot + len(S) <= cmax:
                nxt = cur + [S]
                res.append(nxt)
                rec(i, nxt, tot + len(S))
    rec(0, [], 0)
    return res

def order_bfs(T, M):
    sub = T.subgraph(M)
    out, seen = [], set()
    for comp in nx.connected_components(sub):
        comp = sorted(comp)
        for v in nx.bfs_tree(sub, comp[0]):
            out.append(v); seen.add(v)
    out += sorted(M - seen)
    return out

def run2(d, cmax, mlo, mhi, mod=1, rem=0):
    DELTAS = ['d%d' % i for i in range(d)]
    budget = cmax + 3 * d
    stats = Counter(); survivors = []
    for m in range(mlo, min(mhi, budget - 2) + 1):
        tcount = -1
        for T, k in forests(m, budget):
            tcount += 1
            if tcount % mod != rem: continue
            print('progress m=%d tree=%d' % (m, tcount), file=sys.stderr, flush=True)
            nodes = sorted(T.nodes())
            adjT = {v: set(T[v]) for v in nodes}
            degT = dict(T.degree())
            for Mmask in range(1, 1 << m):
                M = {nodes[i] for i in range(m) if Mmask >> i & 1}
                FixC = set(nodes) - M
                if any(len(adjT[x] & FixC) > 1 for x in M): continue
                Ml = order_bfs(T, M)
                stats['(T,M) pairs'] += 1
                psis = []
                def dfs(i, psi, used):
                    if i == len(Ml):
                        psis.append(dict(psi)); return
                    x = Ml[i]
                    for img in Ml + DELTAS:
                        if img == x or img in used: continue
                        if img in M and (adjT[x] & FixC) != (adjT[img] & FixC):
                            continue
                        ok = True
                        for y in Ml[:i]:
                            iy = psi[y]
                            if img in M and iy in M:
                                if ((iy in adjT[img]) != (y in adjT[x])):
                                    ok = False; break
                        if not ok: continue
                        psi[x] = img; used.add(img)
                        dfs(i + 1, psi, used)
                        del psi[x]; used.discard(img)
                dfs(0, {}, set())
                for psi in psis:
                    stats['psi maps'] += 1
                    Eset = {psi[x] for x in M if psi[x] in DELTAS}
                    # forced structures from pullback through exits
                    FS = defaultdict(set); FA = defaultdict(set)
                    dd_e, dd_a = set(), set()
                    for x in M:
                        e = psi[x]
                        if e not in DELTAS: continue
                        for y in nodes:
                            if y == x: continue
                            a = y in adjT[x]
                            img = psi[y] if y in M else y
                            if img in DELTAS:
                                (dd_e if a else dd_a).add(frozenset((e, img)))
                            else:
                                (FS if a else FA)[img].add(e)
                    if dd_e & dd_a: continue
                    if any(FS[v] & FA[v] for v in list(FS)): continue
                    cands = []
                    for sz in range(1, cmax + 1):
                        for S in combinations(nodes, sz):
                            Sf = frozenset(S)
                            if all((x in Sf) == (psi[x] in Sf)
                                   for x in M if psi[x] not in DELTAS):
                                cands.append(Sf)
                    for ratt in collections_of(cands, cmax):
                        stats['(psi,r)'] += 1
                        rdeg = Counter()
                        for S in ratt:
                            for x in S: rdeg[x] += 1
                        if any(degT[x] + rdeg[x] > 3 for x in nodes): continue
                        s = {x: 3 - degT[x] - rdeg[x] for x in nodes}
                        if any(len(FS[v]) > s[v] for v in nodes): continue
                        rd_e = set()
                        for i, S in enumerate(ratt):
                            for x in M:
                                if psi[x] in DELTAS and x in S:
                                    rd_e.add(('r%d' % i, psi[x]))
                        forced_cnt = Counter()
                        for (_, dl) in rd_e: forced_cnt[dl] += 1
                        for p in dd_e:
                            for dl in p: forced_cnt[dl] += 1
                        if sum(s.values()) + sum(forced_cnt.values()) > 3 * d:
                            continue
                        # base graph and distances for slot pruning
                        adjB = {v: set(adjT[v]) for v in nodes}
                        for i, S in enumerate(ratt):
                            lbl = 'r%d' % i
                            adjB[lbl] = set(S)
                            for x in S: adjB[x].add(lbl)
                        if not girth_ok(adjB): continue
                        Bx = nx.Graph()
                        for v, nb in adjB.items():
                            Bx.add_node(v); Bx.add_edges_from((v, w) for w in nb)
                        dist = dict(nx.all_pairs_shortest_path_length(Bx))
                        verts = sorted(nodes, key=lambda v: -s[v])
                        dsets = {dl: set() for dl in DELTAS}
                        cap = Counter(forced_cnt)
                        pair_used = set()
                        def bt(i):
                            if i == len(verts):
                                finish(); return
                            x = verts[i]
                            base = FS[x]
                            free = [dl for dl in DELTAS
                                    if dl not in base and dl not in FA[x]]
                            need = s[x] - len(base)
                            if need < 0 or need > len(free): return
                            for extra in combinations(free, need):
                                chosen = tuple(base) + extra
                                ok = True; newpairs = []
                                for dl in chosen:
                                    if cap[dl] >= 3: ok = False; break
                                    for y in dsets[dl]:
                                        if dist[x].get(y, 99) < 3: ok = False; break
                                        fp = frozenset((x, y))
                                        if fp in pair_used or fp in newpairs:
                                            ok = False; break
                                        newpairs.append(fp)
                                    if not ok: break
                                if not ok: continue
                                for dl in chosen:
                                    dsets[dl].add(x); cap[dl] += 1
                                for fp in newpairs: pair_used.add(fp)
                                bt(i + 1)
                                for dl in chosen:
                                    dsets[dl].remove(x); cap[dl] -= 1
                                for fp in newpairs: pair_used.discard(fp)
                        def finish():
                            stats['complete assignments'] += 1
                            adjK = {v: set(nb) for v, nb in adjB.items()}
                            for dl in DELTAS:
                                if dsets[dl]:
                                    adjK.setdefault(dl, set())
                                    for x in dsets[dl]:
                                        adjK[dl].add(x); adjK[x].add(dl)
                            for (rl, dl) in rd_e:
                                adjK.setdefault(rl, set()).add(dl)
                                adjK.setdefault(dl, set()).add(rl)
                            for p in dd_e:
                                a_, b_ = tuple(p)
                                adjK.setdefault(a_, set()).add(b_)
                                adjK.setdefault(b_, set()).add(a_)
                            # exit-degree completion (fixpoint)
                            okc = True
                            tgt = [dl for dl in DELTAS if dl not in Eset]
                            changed = True
                            while changed and okc:
                                changed = False
                                for e in sorted(Eset):
                                    known_e = len(adjK.get(e, set()))
                                    r_ = 3 - known_e
                                    if r_ < 0: okc = False; break
                                    free_t = [dl for dl in tgt
                                              if dl not in adjK.get(e, set())]
                                    if r_ > len(free_t): okc = False; break
                                    if r_ == len(free_t) and r_ > 0:
                                        for dl in free_t:
                                            adjK.setdefault(e, set()).add(dl)
                                            adjK.setdefault(dl, set()).add(e)
                                        changed = True
                            if okc:
                                for dl in DELTAS:
                                    if len(adjK.get(dl, set())) > 3: okc = False
                            if not okc:
                                stats['kill exit completion'] += 1; return
                            if not girth_ok(adjK):
                                stats['kill girth'] += 1; return
                            stats['pre-IPR survivors'] += 1
                            if not Eset: stats['pre-IPR (E empty)'] += 1
                            if ipr_verdict(adjK, set(nodes)) is None:
                                stats['kill IPR faces'] += 1
                            else:
                                stats['SURVIVORS'] += 1
                                survivors.append((m, sorted(T.edges()), dict(psi),
                                                  [sorted(S) for S in ratt],
                                                  {dl: sorted(dsets[dl]) for dl in DELTAS}))
                        bt(0)
    print('--- (2prime v2) d=%d c<=%d m=[%d,%d] mod=%d rem=%d ---'
          % (d, cmax, mlo, min(mhi, budget - 2), mod, rem))
    for key in sorted(stats): print('  %-30s %d' % (key, stats[key]))
    if survivors:
        print('  !!! survivors:')
        for sv in survivors: print('   ', sv)
    else:
        print('  no survivors in this range.')
    return survivors

if __name__ == '__main__':
    a = list(map(int, sys.argv[1:]))
    d, cmax, mlo, mhi = a[:4]
    mod, r_ = (a[4], a[5]) if len(a) >= 6 else (1, 0)
    run2(d, cmax, mlo, mhi, mod, r_)
