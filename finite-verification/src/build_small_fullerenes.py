import networkx as nx
from search_fullerene_localised import fullerene_ok, is_ipr, faces_of, localised_pa

def drum(m, bands):
    """m-fold: top m-gon, then `bands` pairs of (A,C) layers of m vertices in the
    dodecahedron zigzag pattern, then bottom m-gon."""
    G=nx.Graph()
    T=[('T',i) for i in range(m)]; B=[('B',i) for i in range(m)]
    for i in range(m):
        G.add_edge(T[i],T[(i+1)%m]); G.add_edge(B[i],B[(i+1)%m])
    layers=[]
    for b in range(bands):
        layers.append([('A%d'%b,i) for i in range(m)])
        layers.append([('C%d'%b,i) for i in range(m)])
    # top m-gon spokes to first layer
    for i in range(m): G.add_edge(T[i], layers[0][i])
    # zigzag between consecutive layers (dodeca pattern): A[i]~C[i],C[i-1]; C[i]~A[i],A[i+1]
    for r in range(len(layers)-1):
        A,C=layers[r],layers[r+1]
        for i in range(m):
            G.add_edge(A[i],C[i])
            if r%2==0: G.add_edge(A[i],C[(i-1)%m])
            else:      G.add_edge(A[i],C[(i+1)%m])
    # last layer spokes to bottom m-gon
    for i in range(m): G.add_edge(layers[-1][i], B[i])
    return G

if __name__=="__main__":
    found=[]
    for m in range(5,8):
        for bands in range(1,6):
            G=drum(m,bands)
            if all(d==3 for _,d in G.degree()) and fullerene_ok(G):
                fs=sorted(set(len(f) for f in faces_of(G)))
                found.append((m,bands,G.number_of_nodes(),is_ipr(G)))
                print("drum m=%d bands=%d -> C%d fullerene, IPR=%s faces=%s"%(
                    m,bands,G.number_of_nodes(),is_ipr(G),fs))
    print("\nvalid fullerenes found:",len(found))
