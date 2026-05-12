import Pkg
Pkg.activate(".")

using GraphIO.Graph6
using Graphs
using NautyGraphs
using Combinatorics
using DataFrames
using CSV
using Base.Threads
using SHA

# Function to generate 128-bit hash using SHA256
function hash_sha128(x...)
    io = IOBuffer()
    for val in x
        if val isa AbstractArray
            for item in val
                write(io, htol(item))
            end
        else
            write(io, htol(val))
        end
    end
    return reinterpret(UInt128, sha256(take!(io)))[1]
end

# Get grpsize and 128-bit hash of canonical form
function canonical_info(g::DenseNautyGraph)
    canong, _, _, statistics = NautyGraphs._nauty(g)
    grpsize = Float64(statistics.grpsize1) * 10.0^statistics.grpsize2
    h = hash_sha128(canong.words)
    return grpsize, h
end

@inline function check_duplicates(arr::AbstractVector{UInt128})
    @inbounds for i in 1:length(arr)-1
        if arr[i] == arr[i + 1]
            return true
        end
    end
    return false
end

function asymmetric_depth(g6_string::String; max_depth=5)
    g_simple = Graph6._g6StringToGraph(g6_string)
    A = Graphs.adjacency_matrix(g_simple)
    g = NautyGraph(A)
    n = nv(g)
    
    # Check if the graph itself is symmetric
    grpsize, _ = canonical_info(g)
    if grpsize > 1.0
        return 0
    end

    for i in 1:max_depth
        canonical_table = UInt128[]
        
        for comb in combinations(1:n, n - i)
            subg_simple, _ = Graphs.induced_subgraph(g_simple, comb)
            A_sub = Graphs.adjacency_matrix(subg_simple)
            g_sub = NautyGraph(A_sub)
            
            grpsize_sub, h = canonical_info(g_sub)
            if grpsize_sub > 1.0
                return i
            end
            push!(canonical_table, h)
        end
        
        sort!(canonical_table)
        if check_duplicates(canonical_table)
            return i
        end
    end
    return max_depth + 1
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
    for n_c in 60:2:100
        g6_c_list, g6_t_list = get_ipr_fullerenes(n_c)
        n = min(length(g6_c_list), length(g6_t_list))
        for i in 1:n
            push!(all_tasks, (n_c=n_c, gc=g6_c_list[i], gt=g6_t_list[i]))
        end
    end
    
    total = length(all_tasks)
    println("Total IPR fullerenes to process: $total")
    
    results = Vector{Any}(nothing, total)
    
    # Depth limit set to 4
    Threads.@threads for i in 1:total
        task = all_tasks[i]
        dc = asymmetric_depth(task.gc, max_depth=4)
        dt = asymmetric_depth(task.gt, max_depth=4)
        results[i] = (n_c=task.n_c, gc=task.gc, gt=task.gt, depth_cubic=dc, depth_tri=dt)
        if i % 50 == 0
            println("Completed $i / $total")
        end
    end
    
    df = DataFrame(n_c = Int[], fullerene_g6 = String[], dual_g6 = String[], depth_fullerene = Int[], depth_dual = Int[])
    for res in results
        if res !== nothing
            push!(df, (n_c=res.n_c, fullerene_g6=res.gc, dual_g6=res.gt, depth_fullerene=res.depth_cubic, depth_dual=res.depth_tri))
        end
    end
    
    if !isdir("results")
        mkdir("results")
    end
    
    CSV.write("results/ipr_fullerenes_depths_128bit.csv", df)
    println("Results saved to results/ipr_fullerenes_depths_128bit.csv")
end

main()