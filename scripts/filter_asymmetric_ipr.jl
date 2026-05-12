using GraphIO.Graph6
using Graphs
using NautyGraphs

function is_asymmetric(g6_string::AbstractString)
    g6 = Graph6._g6StringToGraph(g6_string)
    A = Graphs.adjacency_matrix(g6)
    g = NautyGraph(A)
    grpsize::Int64, _, _ = NautyGraphs.nauty(g, false)
    return grpsize == 1
end

function main()
    for line in eachline(stdin)
        s = strip(line)
        if isempty(s) || s == ">>graph6<<"
            continue
        end
        if is_asymmetric(s)
            println(s)
        end
    end
end

main()
