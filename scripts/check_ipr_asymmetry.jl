
import Pkg
Pkg.activate(".")

using GraphIO.Graph6
using Graphs
using NautyGraphs

function is_asymmetric(g6_string::AbstractString)
    try
        g6 = Graph6._g6StringToGraph(String(g6_string))
        A = Graphs.adjacency_matrix(g6)
        g = NautyGraph(A)
        grpsize::Int64, _, _ = NautyGraphs.nauty(g, false)
        return grpsize == 1
    catch
        return false
    end
end

function main()
    println("n_c,ipr_count,asymmetric_count")
    # Store asymmetric ones to a file for later
    open("asymmetric_ipr_fullerenes.txt", "w") do asym_io
        for n in 60:2:100
            # Call buckygen to get IPR fullerenes
            # Use -q to suppress stderr if possible, or just ignore it in the script
            cmd = `./buckygen-1.1/buckygen-1.1/buckygen -I -d -g $(n)d`
            
            ipr_count = 0
            asymmetric_count = 0
            
            # We want to read stdout only. 
            # In Julia, open(cmd, "r") reads stdout.
            open(cmd, "r") do io
                for line in eachline(io)
                    s = strip(line)
                    if isempty(s) || s == ">>graph6<<"
                        continue
                    end
                    # Graph6 strings for these sizes don't start with 'N' or '1 '
                    if length(s) > 0 && !startswith(s, "Number of") && !startswith(s, "1 fullerenes written")
                        ipr_count += 1
                        if is_asymmetric(s)
                            asymmetric_count += 1
                            println(asym_io, "$n $s")
                        end
                    end
                end
            end
            if ipr_count > 0
                println("$n,$ipr_count,$asymmetric_count")
            end
        end
    end
end

main()
