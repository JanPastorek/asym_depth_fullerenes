import itertools, networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher

def faces_of(G):
    pl,emb=nx.check_planarity(G); seen=set(); fs=[]
    for u in emb:
        for v in emb[u]:
            if (u,v) in seen: continue
            fs.append(emb.traverse_face(u,v,mark_half_edges=seen))
    return fs
def fullerene_ok(G):
    if any(d!=3 for _,d in G.degree()): return False
    pl,_=nx.check_planarity(G)
    if not pl or nx.node_connectivity(G)<3 or nx.girth(G)<5: return False
    return set(len(f) for f in faces_of(G))<={5,6}
def is_ipr(G):
    pe=[{frozenset((f[i],f[(i+1)%5])) for i in range(5)} for f in faces_of(G) if len(f)==5]
    return not any(a&b for a,b in itertools.combinations(pe,2))

def uf_components(nodes, edges):
    p={x:x for x in nodes}
    def find(x):
        while p[x]!=x: p[x]=p[p[x]]; x=p[x]
        return x
    for a,b in edges:
        ra,rb=find(a),find(b)
        if ra!=rb: p[ra]=rb
    comp={}
    for x in nodes: comp.setdefault(find(x),set()).add(x)
    return list(comp.values())

def localised_pa(G, ks=(3,2,1), cmax_C=9):
    V=list(G.nodes())
    adj={v:set(G[v]) for v in V}
    Eall=[tuple(sorted(e)) for e in G.edges()]
    for k in ks:
        cbound=5-k
        for SD in itertools.combinations(V,k):
            SDs=set(SD)
            Hn=[v for v in V if v not in SDs]
            E=[e for e in Eall if e[0] not in SDs and e[1] not in SDs]
            for r in range(1,cbound+1):
                for X in itertools.combinations(range(len(E)),r):
                    Xset=set(X)
                    keep=[E[i] for i in range(len(E)) if i not in Xset]
                    comps=uf_components(Hn,keep)
                    if len(comps)<2: continue
                    for C in comps:
                        if not (2<=len(C)<=cmax_C): continue
                        R=set(Hn)-C
                        if len(R)<6: continue
                        cut=sum(1 for a,b in E if (a in C)!=(b in C))
                        if cut>cbound: continue
                        GC=G.subgraph(C)
                        NR={c:frozenset(adj[c]&R) for c in C}
                        for mp in GraphMatcher(GC,GC).isomorphisms_iter():
                            if all(mp[c]==c for c in C): continue
                            if all(NR[c]==NR[mp[c]] for c in C):
                                return dict(k=k,SD=SD,C=tuple(sorted(C)),Rsize=len(R),
                                            interface=cut,alpha={a:mp[a] for a in C if mp[a]!=a})
    return None

if __name__=="__main__":
    import sys
    name=sys.argv[1] if len(sys.argv)>1 else "dodecahedral"
    G=getattr(nx,name+"_graph")() if name in ("dodecahedral",) else None
    print(name,"fullerene:",fullerene_ok(G),"IPR:",is_ipr(G))
    import time; t=time.time()
    print("witness:",localised_pa(G, ks=(3,)))
    print("k=3 time %.1fs"%(time.time()-t))
