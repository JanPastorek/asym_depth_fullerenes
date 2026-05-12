import Pkg
Pkg.activate(".")

using GraphIO.Graph6
using Graphs
using NautyGraphs
using Combinatorics
using DataFrames
using CSV
using Base.Threads

function asymmetric_depth(g::NautyGraph)
    n = nv(g)
    
    # Check if the graph itself is symmetric
    _, autg = nauty(g)
    if autg.n > 1.0
        return 0
    end

    # i is the number of vertices to remove
    for i in 1:n-1
        # Use a vector for hashes to be memory efficient and allow sorting
        canonical_table = UInt128[]
        
        for comb in combinations(1:n, n - i)
            # NautyGraphs.induced_subgraph returns a NautyGraph directly
            subg, _ = Graphs.induced_subgraph(g, comb)
            
            # Check for symmetry in the subgraph
            _, autg_sub = nauty(subg)
            if autg_sub.n > 1.0
                return i
            end
            
            # Get 128-bit canonical hash
            h = canonical_id(subg)
            push!(canonical_table, h)
        end
        
        # Check for duplicates in the canonical table (meaning isomorphic subgraphs)
        sort!(canonical_table)
        for j in 1:length(canonical_table)-1
            if canonical_table[j] == canonical_table[j+1]
                return i
            end
        end
    end
    return n - 1
end

function is_graph6(s::AbstractString)
    return length(s) > 0 && (63 <= Int(s[1]) <= 126) && !occursin(" ", s) && !occursin("Number", s)
end

function get_ipr_fullerenes(n_c)
    cmd_t = `./buckygen-1.1/buckygen-1.1/buckygen -I -g $(n_c)d`
    lines_t = try readlines(cmd_t) catch; String[] end
    g6_t = String[String(strip(l)) for l in lines_t if is_graph6(strip(l))]
    
    cmd_c = `./buckygen-1.1/buckygen-1.1/buckygen -I -d -g $(n_c)d`
    lines_c = try readlines(cmd_c) catch; String[] end
    g6_c = String[String(strip(l)) for l in lines_c if is_graph6(strip(l))]
    
    return g6_c, g6_t
end

function main()
    all_tasks = []
    println("Collecting IPR fullerenes...")
    flush(stdout)
    for n_c in 60:2:100
        g6_c_list, g6_t_list = get_ipr_fullerenes(n_c)
        n = min(length(g6_c_list), length(g6_t_list))
        for i in 1:n
            push!(all_tasks, (n_c=n_c, gc=g6_c_list[i], gt=g6_t_list[i]))
        end
    end
    
    total = length(all_tasks)
    println("Total IPR fullerenes to process: $total")
    flush(stdout)
    
    output_path = "results/ipr_fullerenes_depths_v3.csv"
    if !isdir("results")
        mkdir("results")
    end
    
    results_df = if isfile(output_path)
        CSV.read(output_path, DataFrame)
    else
        DataFrame(n_c = Int[], fullerene_g6 = String[], dual_g6 = String[], depth_fullerene = Int[], depth_dual = Int[])
    end
    
    processed_g6 = Set(results_df.fullerene_g6)
    df_lock = ReentrantLock()
    
    println("Starting depth calculations with $(nthreads()) threads...")
    flush(stdout)
    
    @threads for i in 1:total
        task = all_tasks[i]
        
        if task.gc in processed_g6
            continue
        end
        
        # Cubic asymmetric depth
        gc_simple = Graph6._g6StringToGraph(task.gc)
        g_cubic = NautyGraph(Graphs.adjacency_matrix(gc_simple))
        dc = asymmetric_depth(g_cubic)
        
        # Triangulation asymmetric depth
        gt_simple = Graph6._g6StringToGraph(task.gt)
        g_tri = NautyGraph(Graphs.adjacency_matrix(gt_simple))
        dt = asymmetric_depth(g_tri)
        
        lock(df_lock) do
            push!(results_df, (n_c=task.n_c, fullerene_g6=task.gc, dual_g6=task.gt, depth_fullerene=dc, depth_dual=dt))
            if size(results_df, 1) % 10 == 0
                CSV.write(output_path, results_df)
                println("Progress: $(size(results_df, 1)) / $total")
                flush(stdout)
            end
        end
    end
    
    CSV.write(output_path, results_df)
    println("Done! Results saved to $output_path")
    flush(stdout)
end

main()