import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV
using DataFrames
using Glob

function main()
    if length(ARGS) < 2
        println("Usage: julia scripts/merge_asym_depth_results.jl <input_glob> <output_csv>")
        exit(1)
    end

    input_glob = ARGS[1]
    output_csv = ARGS[2]

    files = sort(glob(input_glob))
    if isempty(files)
        error("No files matched glob: $input_glob")
    end

    frames = DataFrame[]
    for f in files
        push!(frames, CSV.read(f, DataFrame))
    end

    df = vcat(frames...)
    unique!(df, [:n_c, :fullerene_g6, :dual_g6])
    sort!(df, [:n_c, :depth_fullerene, :depth_dual])

    outdir = dirname(output_csv)
    if !isempty(outdir) && !isdir(outdir)
        mkpath(outdir)
    end
    CSV.write(output_csv, df)

    println("Merged files: $(length(files))")
    println("Rows: $(nrow(df))")
    println("Saved: $output_csv")
end

main()