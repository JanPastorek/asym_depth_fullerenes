import Pkg
Pkg.activate(".")

using GraphIO.Graph6
using Graphs
using NautyGraphs
using Combinatorics
using DataFrames
using CSV
using Base.Threads

function asymmetric_depth(g::NautyGraph; max_depth=3)
    n = g.n_vertices
    grpsize::Int64, _, _ = NautyGraphs.nauty(g, false)
    if grpsize > 1
        return 0
    end

    for i in 1:max_depth
        canonical_table = Set{UInt64}()
        for comb in combinations(1:n, n - i)
            subg, _ = Graphs.induced_subgraph(g, comb)
            g_sub = NautyGraph(subg)
            
            grpsize_sub, _, _, canon_hash_sub::UInt64 = NautyGraphs._nautyhash(g_sub)
            if grpsize_sub > 1 || (canon_hash_sub in canonical_table)
                return i
            end
            push!(canonical_table, canon_hash_sub)
        end
    end
    return max_depth + 1 # Means > max_depth
end

function is_graph6(s::AbstractString)
    return length(s) > 0 && (63 <= Int(s[1]) <= 126) && !occursin(" ", s) && !occursin("Number", s)
end

function get_ipr_fullerenes(n_c)
    cmd_t = `./buckygen-1.1/buckygen-1.1/buckygen -I -g $(n_c)d`
    lines_t = readlines(cmd_t)
    g6_t = String[String(strip(l)) for l in lines_t if is_graph6(strip(l))]
    
    cmd_c = `./buckygen-1.1/buckygen-1.1/buckygen -I -d -g $(n_c)d`
    lines_c = readlines(cmd_c)
    g6_c = String[String(strip(l)) for l in lines_c if is_graph6(strip(l))]
    
    return g6_c, g6_t
end

function main()
    all_tasks = []
    for n_c in 20:2:100
        g6_c_list, g6_t_list = get_ipr_fullerenes(n_c)
        # Ensure we have the same number of graphs
        n = min(length(g6_c_list), length(g6_t_list))
        for i in 1:n
            push!(all_tasks, (n_c, g6_c_list[i], g6_t_list[i]))
        end
    end
    
    total = length(all_tasks)
    println("Total fullerenes to process: $total")
    
    results = Vector{Any}(nothing, total)
    
    # Process tasks. Using threads if available.
    try
        Threads.@threads for i in 1:total
            n_c, gc, gt = all_tasks[i]
            
            # Cubic depth
            Ac = Graphs.adjacency_matrix(Graph6._g6StringToGraph(gc))
            g_cubic = NautyGraph(Ac)
            dc = asymmetric_depth(g_cubic, max_depth=3)
            
            # Triangulation depth
            At = Graphs.adjacency_matrix(Graph6._g6StringToGraph(gt))
            g_tri = NautyGraph(At)
            dt = asymmetric_depth(g_tri, max_depth=3)
            
            results[i] = (n_c, gc, gt, dc, dt)
            if i % 50 == 0
                println("Completed $i / $total")
            end
        end
    catch e
        println("Parallel execution failed, falling back to sequential: $e")
        for i in 1:total
            if results[i] === nothing
                n_c, gc, gt = all_tasks[i]
                Ac = Graphs.adjacency_matrix(Graph6._g6StringToGraph(gc))
                g_cubic = NautyGraph(Ac)
                dc = asymmetric_depth(g_cubic, max_depth=3)
                At = Graphs.adjacency_matrix(Graph6._g6StringToGraph(gt))
                g_tri = NautyGraph(At)
                dt = asymmetric_depth(g_tri, max_depth=3)
                results[i] = (n_c, gc, gt, dc, dt)
            end
        end
    end
    
    df = DataFrame(n_c = Int[], fullerene_g6 = String[], dual_g6 = String[], fullerene_depth = Int[], dual_depth = Int[])
    for res in results
        if res !== nothing
            push!(df, (n_c=res[1], fullerene_g6=res[2], dual_g6=res[3], fullerene_depth=res[4], dual_depth=res[5]))
        end
    end
    
    if !isdir("results")
        mkdir("results")
    end
    
    CSV.write("results/ipr_depths.csv", df)
    println("Done! Results saved to results/ipr_depths.csv")
end

main()