import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using GraphIO.Graph6
using Graphs
using NautyGraphs
using Combinatorics
using DataFrames
using CSV
using Base.Threads
using SHA

# Compute a stable 128-bit hash from canonical nauty words.
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

function asymmetric_depth(g6_string::String; max_depth::Int=4)
    g_simple = Graph6._g6StringToGraph(g6_string)
    A = Graphs.adjacency_matrix(g_simple)
    g = NautyGraph(A)
    n = nv(g)

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
    return !isempty(s) && (63 <= Int(s[1]) <= 126) && !occursin(" ", s) && !occursin("Number", s)
end

function get_ipr_fullerenes(n_c::Int)
    candidates = String[
        get(ENV, "BUCKYGEN_BIN", ""),
        joinpath(@__DIR__, "..", "buckygen-1.1", "buckygen"),
        joinpath(@__DIR__, "..", "..", "buckygen-1.1", "buckygen")
    ]
    base = ""
    for cand in candidates
        if !isempty(cand) && isfile(cand)
            base = cand
            break
        end
    end

    if isempty(base)
        error("buckygen executable not found. Set BUCKYGEN_BIN or place it at project/asym_depth/buckygen-1.1/buckygen")
    end

    cmd_t = `$base -I -g $(n_c)d`
    cmd_c = `$base -I -d -g $(n_c)d`

    lines_t = try
        readlines(cmd_t)
    catch err
        println("n_c=$n_c: buckygen failed for dual graphs: $err")
        String[]
    end
    lines_c = try
        readlines(cmd_c)
    catch err
        println("n_c=$n_c: buckygen failed for cubic graphs: $err")
        String[]
    end

    g6_t = String[String(strip(l)) for l in lines_t if is_graph6(strip(l))]
    g6_c = String[String(strip(l)) for l in lines_c if is_graph6(strip(l))]
    return g6_c, g6_t
end

function parse_int_arg(s::AbstractString, name::AbstractString)
    try
        return parse(Int, s)
    catch
        error("Invalid $name: $s")
    end
end

function main()
    if length(ARGS) < 3
        println("Usage: julia scripts/calculate_fullerene_asym_depth_hpc.jl <n_start> <n_end> <output_csv> [max_depth]")
        exit(1)
    end

    n_start = parse_int_arg(ARGS[1], "n_start")
    n_end = parse_int_arg(ARGS[2], "n_end")
    output_csv = ARGS[3]
    max_depth = length(ARGS) >= 4 ? parse_int_arg(ARGS[4], "max_depth") : 4

    if n_start > n_end
        error("n_start must be <= n_end")
    end

    n_values = [n for n in n_start:2:n_end if iseven(n)]
    if isempty(n_values)
        error("No even n_c values in range [$n_start, $n_end]")
    end

    all_tasks = Vector{NamedTuple{(:n_c, :gc, :gt), Tuple{Int, String, String}}}()
    for n_c in n_values
        g6_c_list, g6_t_list = get_ipr_fullerenes(n_c)
        n = min(length(g6_c_list), length(g6_t_list))
        if n == 0
            println("n_c=$n_c: no IPR fullerenes found")
            continue
        end
        if length(g6_c_list) != length(g6_t_list)
            println("n_c=$n_c: warning mismatch cubic=$(length(g6_c_list)) dual=$(length(g6_t_list)); using first $n pairs")
        end
        for i in 1:n
            push!(all_tasks, (n_c=n_c, gc=g6_c_list[i], gt=g6_t_list[i]))
        end
        println("n_c=$n_c: tasks=$n")
    end

    total = length(all_tasks)
    println("Total tasks: $total")
    println("Threads: $(nthreads())")

    results = Vector{Any}(nothing, total)
    @threads for i in 1:total
        task = all_tasks[i]
        dc = asymmetric_depth(task.gc, max_depth=max_depth)
        dt = asymmetric_depth(task.gt, max_depth=max_depth)
        results[i] = (n_c=task.n_c, fullerene_g6=task.gc, dual_g6=task.gt, depth_fullerene=dc, depth_dual=dt)
        if i % 25 == 0
            println("Completed $i / $total")
        end
    end

    df = DataFrame(n_c=Int[], fullerene_g6=String[], dual_g6=String[], depth_fullerene=Int[], depth_dual=Int[])
    for r in results
        if r !== nothing
            push!(df, r)
        end
    end

    outdir = dirname(output_csv)
    if !isempty(outdir) && !isdir(outdir)
        mkpath(outdir)
    end
    CSV.write(output_csv, df)
    println("Saved: $output_csv")
end

main()