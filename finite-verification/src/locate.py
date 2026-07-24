import sys, subprocess
from pcode import parse_planarcode
from search_fullerene_localised import is_ipr, faces_of

def sep_count(G):
    faces={frozenset(f) for f in faces_of(G) if len(f)==5}
    A={v:set(G[v]) for v in G}; nd=sorted(G.nodes()); idx={v:i for i,v in enumerate(nd)}
    allf=set()
    for s in nd:
        st=[(s,(s,))]
        while st:
            v,p=st.pop()
            for w in A[v]:
                if w==s and len(p)==5:
                    if idx[p[1]]<idx[p[-1]]: allf.add(frozenset(p))
                elif idx[w]>idx[s] and w not in p and len(p)<5: st.append((w,p+(w,)))
    return len(allf-faces)

for V in [int(x) for x in sys.argv[1:]]:
    n=(V+4)//2
    gs=parse_planarcode(subprocess.run(["/tmp/buckygen/buckygen","-d",str(n)],capture_output=True).stdout)
    sep=sum(1 for G in gs if sep_count(G)>0)
    print("C%d: %d isomers, %d with >=1 separating pentagon"%(V,len(gs),sep),flush=True)
