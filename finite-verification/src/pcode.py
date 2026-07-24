import sys, networkx as nx
def parse_planarcode(data):
    # skip header >>planar_code<< ... up to and including the '<<'
    hdr=b'>>planar_code'
    i=data.find(hdr)
    # header ends at first '<<' after start
    j=data.find(b'<<', i+len(hdr))
    p=j+2
    graphs=[]
    n=len(data)
    while p<n:
        order=data[p]; p+=1
        if order==0: break
        G=nx.Graph(); G.add_nodes_from(range(1,order+1))
        ok=True
        for v in range(1,order+1):
            while p<n and data[p]!=0:
                G.add_edge(v,data[p]); p+=1
            if p>=n: ok=False; break
            p+=1  # skip 0
        graphs.append(G)
    return graphs

if __name__=="__main__":
    data=sys.stdin.buffer.read()
    gs=parse_planarcode(data)
    from search_fullerene_localised import fullerene_ok, is_ipr
    nfull=sum(1 for G in gs if fullerene_ok(G))
    print("parsed graphs:",len(gs),"| valid fullerenes:",nfull,
          "| sizes:",sorted(set(G.number_of_nodes() for G in gs)))
    if gs:
        G=gs[0]; print("first: n=%d fullerene=%s IPR=%s"%(G.number_of_nodes(),fullerene_ok(G),is_ipr(G)))
