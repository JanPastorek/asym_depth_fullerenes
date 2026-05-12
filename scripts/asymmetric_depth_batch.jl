using GraphIO.Graph6
using Graphs
using NautyGraphs
using Combinatorics

function asymmetric_depth(g6_string::String)
    g6 = Graph6._g6StringToGraph(g6_string)
    A = Graphs.adjacency_matrix(g6)
    g = NautyGraph(A)
    canonical_table = Set{UInt64}()
    n = g.n_vertices
    grpsize::Int64, _, _ = NautyGraphs.nauty(g, false)
    if grpsize > 1
        return 0
    end

    @inbounds for i in 1:3
        @inbounds for comb in combinations(1:n, n - i)
            subg, _ = Graphs.induced_subgraph(g, comb)
            grpsize, _, _, canon_hash_sub::UInt64 = NautyGraphs._nautyhash(subg)
            if grpsize > 1 || (canon_hash_sub in canonical_table)
                return i
            end
            push!(canonical_table, canon_hash_sub)
        end
        empty!(canonical_table)
    end
    return 4
end

function main()
    if length(ARGS) < 2
        println("Usage: julia asymmetric_depth_batch.jl <input_g6> <output_csv>")
        exit(1)
    end
    input_path = ARGS[1]
    output_path = ARGS[2]

    open(output_path, "w") do io
        println(io, "g6,asymmetric_depth")
        for line in eachline(input_path)
            g6 = String(strip(line))
            if isempty(g6)
                continue
            end
            depth = asymmetric_depth(g6)
            println(io, g6, ",", depth)
        end
    end
end

main()
