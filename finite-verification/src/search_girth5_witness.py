import random, networkx as nx, sys
from itertools import combinations
from networkx.algorithms.isomorphism import GraphMatcher
C = ['x','u','w','v','y']
FRAG = [('x','u'),('u','w'),('w','v'),('v','y'),('rho','x'),('rho','y'),
        ('s1','u'),('s1','y'),('s2','v'),('s2','x'),('s3','w')]
PSI = {'x':'y','y':'x','u':'v','v':'u','w':'w'}
DEL = {'s1','s2','s3'}
def is_pa(G, phi):
    for a,b in combinations(sorted(phi),2):
        if G.has_edge(a,b) != G.has_edge(phi[a],phi[b]): return False
    return True
def naut(G, cap=2):
    n=0
    for _ in GraphMatcher(G,G).isomorphisms_iter():
        n+=1
        if n>=cap: return n
    return n
def attempt(m, rng):
    extra=['e%d'%i for i in range(m)]
    G=nx.Graph(); G.add_edges_from(FRAG); G.add_nodes_from(extra)
    stubs=['rho','s1','s2','s3','s3']+[v for v in extra for _ in range(3)]
    rng.shuffle(stubs)
    if len(stubs)%2: return None
    for i in range(0,len(stubs),2):
        a,b=stubs[i],stubs[i+1]
        if a==b or G.has_edge(a,b): return None
        G.add_edge(a,b)
    if any(d!=3 for _,d in G.degree()): return None
    if nx.girth(G)<5: return None
    if not nx.is_connected(G) or nx.node_connectivity(G)<3: return None
    phi={v:PSI[v] for v in C}
    for v in G.nodes():
        if v not in C and v not in DEL: phi[v]=v
    if not is_pa(G,phi): return None
    return G,phi
rng=random.Random(int(sys.argv[1]) if len(sys.argv)>1 else 1)
for m in range(5,30,2):
    for _ in range(3000):
        r=attempt(m,rng)
        if not r: continue
        G,phi=r
        if naut(G)==1:
            n=G.number_of_nodes()
            print("ASYMMETRIC WITNESS n =",n,"girth",nx.girth(G),
                  "conn",nx.node_connectivity(G))
            print(" k =",n-len(phi),"supp =",sorted(v for v in phi if phi[v]!=v))
            print(" interface =",sum(1 for p,q in G.edges()
                  if (p in C)!=(q in C) and p not in DEL and q not in DEL))
            print(" edges =",sorted(tuple(sorted(e)) for e in G.edges()))
            sys.exit()
    print("m =",m,"none",flush=True)
