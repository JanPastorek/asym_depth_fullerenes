import Pkg
Pkg.activate(".")

using GraphIO.Graph6
using Graphs
using NautyGraphs
using Combinatorics
using TimerOutputs

function asymmetric_depth(g::NautyGraph)
    canonical_table = Set{UInt64}()
    n = g.n_vertices
    grpsize::Int64, _, _ = NautyGraphs.nauty(g, false)
    if grpsize > 1
        return 0
    end

    for i in 1:n-1
        println("Checking depth $i (combinations: $(binomial(n, i)))")
        for comb in combinations(1:n, n - i)
            subg, _ = Graphs.induced_subgraph(g, comb)
            A_sub = Graphs.adjacency_matrix(subg)
            g_sub = NautyGraph(A_sub)
            
            grpsize_sub, _, _, canon_hash_sub::UInt64 = NautyGraphs._nautyhash(g_sub)
            if grpsize_sub > 1 || (canon_hash_sub in canonical_table)
                return i
            end
            push!(canonical_table, canon_hash_sub)
        end
        empty!(canonical_table)
    end
    return n - 1
end

function main()
    # Read the output of buckygen fully to avoid broken pipe
    cmd = `./buckygen-1.1/buckygen-1.1/buckygen -I -d -g 84d`
    lines = readlines(cmd)
    
    g6_cubic = ""
    for line in lines
        s = strip(line)
        if isempty(s) || s == ">>graph6<<" || startswith(s, "Number of") || startswith(s, "1 fullerenes written")
            continue
        end
        # Check if asymmetric
        A = Graphs.adjacency_matrix(Graph6._g6StringToGraph(String(s)))
        g = NautyGraph(A)
        grpsize, _, _ = NautyGraphs.nauty(g, false)
        if grpsize == 1
            g6_cubic = String(s)
            break
        end
    end
    
    if isempty(g6_cubic)
        println("No asymmetric IPR found for n_c=84")
        return
    end
    
    println("Testing n_c=84 fullerene: $g6_cubic")
    A = Graphs.adjacency_matrix(Graph6._g6StringToGraph(g6_cubic))
    g = NautyGraph(A)
    
    to = TimerOutput()
    @timeit to "asymmetric_depth" begin
        depth = asymmetric_depth(g)
        println("Depth: $depth")
    end
    show(to)
    println()
end

main()