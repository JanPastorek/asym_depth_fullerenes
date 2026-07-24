import sys, subprocess
from pcode import parse_planarcode
from search_fullerene_localised import faces_of
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
V=int(sys.argv[1]); mod=int(sys.argv[2]) if len(sys.argv)>2 else 1
n=(V+4)//2
tot=0; sep=0
for res in range(mod):
    arg=[f"./buckygen","-d",str(n)]+([f"{res}/{mod}"] if mod>1 else [])
    data=subprocess.run(arg,capture_output=True).stdout
    for G in parse_planarcode(data):
        tot+=1
        if sep_count(G)>0: sep+=1
    break  # one residue class = a sample of size ~total/mod
print("C%d: sampled %d isomers (1/%d), with separating pentagon: %d"%(V,tot,mod,sep))
