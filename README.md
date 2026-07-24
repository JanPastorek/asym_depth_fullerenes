# Asymmetric depth in fullerenes — code & data

Computational material accompanying the paper

> **Extremal Asymmetric Depth on Surfaces and Hidden Almost-Global Symmetries in IPR Fullerenes**
> Ján Pastorek (Comenius University, Bratislava).

The repository has two independent parts, matching the two kinds of
computational claim in the paper.

## A. Asymmetric-depth census (Julia + HPC)

Computes the asymmetric depth of every IPR fullerene (and its dual) over the
tested range and produces the paper's main table. Built on `NautyGraphs`
(canonical hashing of induced subgraphs); the core routine is
`asymmetric_depth` in `scripts/calculate_ipr_fullerenes_asym_depth_v3.jl`.

| Path | Role |
|---|---|
| `scripts/*.jl` | depth computation, IPR/asymmetry filtering, batching, merging (`calculate_ipr_fullerenes_asym_depth*.jl`, `check_ipr_asymmetry.jl`, `filter_asymmetric_ipr.jl`, `merge_asym_depth_results.jl`, …) |
| `scripts/*.py` | result aggregation and metrics (`aggregate_results.py`, `compute_metrics.py`, `split_g6.py`) |
| `slurm/submit_asym_depth_pipeline.sh` | SLURM array driver for the HPC run (`N_MIN`…`N_MAX`) |
| `buckygen-1.1/` | bundled fullerene generator + `generate_table.py` |
| `results/` | computed depth tables, incl. `ipr_fullerenes_depths_to120.csv` and per-`n` parts under `results/asym_depth_parts/` |
| `Project.toml`, `Manifest.toml` | pinned Julia environment (CSV, Combinatorics, DataFrames, Glob, GraphIO, Graphs, NautyGraphs) |

Reproduce (single machine):

```bash
julia --project=. -e 'import Pkg; Pkg.instantiate()'
julia --project=. scripts/calculate_ipr_fullerenes_asym_depth_v3.jl   # see script header for args
julia --project=. scripts/merge_asym_depth_results.jl 'results/asym_depth_parts/ipr_depth_n*.csv' results/merged.csv
```
For the full range use the SLURM driver in `slurm/`.

**Paper item:** the census table (counts of IPR / asymmetric IPR fullerenes and
the depth-2/3/4 breakdown up to n = 118), and the 47-vertex extremal duals.

## B. Finite verification of the structural theorems (Python)

Self-contained checks (NetworkX only) behind the rigidity results. See
[`finite-verification/README.md`](finite-verification/README.md) for details.

| Path | Paper item |
|---|---|
| `finite-verification/src/search_exit2.py` (engine `search_general.py`) | **Lemma (Computational)** — the acyclic case of the Unified Localisation Theorem; confirms *no configuration record survives* for k ∈ {1,2,3} |
| `finite-verification/src/gen_forests_nauty.py`, `search_forest_case.py` | independent re-run via nauty `gentreeg` |
| `finite-verification/src/verify_witness_paper.py` | **Proposition (girth-5 witness)** — verifies the explicit 20-vertex graph and all its claimed properties |
| `finite-verification/src/{witness,search_girth5_witness}.py` | heuristic witness search |
| `finite-verification/src/{locate,sample,pcode,search_fullerene_localised,build_small_fullerenes}.py` | separating-pentagon census over `buckygen` output; fullerene predicates |

Reproduce everything self-contained:

```bash
cd finite-verification
pip install -r requirements.txt
./run_all.sh          # ~2 min on one core; prints a pass/fail summary
```

## Citing

See [`CITATION.cff`](CITATION.cff). Licensed under the MIT License
([`LICENSE`](LICENSE)).
