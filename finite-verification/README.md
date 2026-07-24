# Finite verification (search & witness)

Computer-assisted proofs and enumerations accompanying the paper

> **Extremal Asymmetric Depth on Surfaces and Hidden Almost-Global Symmetries in IPR Fullerenes**
> Ján Pastorek (Comenius University, Bratislava).

The paper's structural results are partly computer-assisted. This directory
contains the exact programs referenced in the paper's *Computational
methodology* section, so that every finite verification can be reproduced
independently.

All scripts are plain Python 3 and depend only on [NetworkX](https://networkx.org)
(`pip install -r requirements.txt`). Two enumerations additionally call external,
independently trusted binaries from the [`nauty`](https://pallini.di.uniroma1.it/)
suite (`gentreeg`) and [`buckygen`](https://caagt.ugent.be/buckygen/); see
*Prerequisites* below.

## What is verified here

| Paper item | Script(s) | What it does |
|---|---|---|
| **Lemma (Computational), acyclic case** of the Unified Localisation Theorem | `src/search_exit2.py` (engine: `src/search_general.py`) | Exhaustively enumerates the forest-carrier configuration records for deficiency `k ∈ {1,2,3}` and interface budget `c ≤ 5−k`, applies the filters (A1)–(A5), and confirms **no configuration record survives**. |
| Independent re-run of the same case on a trusted tree generator | `src/gen_forests_nauty.py`, `src/search_forest_case.py` | Rebuilds the subcubic carrier trees (`|C| ≤ 9`, `Δ ≤ 3`) with `nauty`'s `gentreeg` and reproduces the zero-survivor result. |
| **Proposition (girth-5 witness)** | `src/verify_witness_paper.py` (exact graph), `src/witness.py`, `src/search_girth5_witness.py` (search) | `verify_witness_paper.py` builds the **explicit 20-vertex graph listed in the paper** and checks every claimed property (cubic, 3-connected, girth 5, trivial Aut, non-planar, rank-17 partial automorphism, support `{x,u,v,y}`, 2-edge interface). The `*witness*` search scripts find such completions heuristically. |
| Absence of separating pentagons / localised configs in non-IPR fullerenes | `src/locate.py`, `src/sample.py`, `src/pcode.py` | Streams `buckygen` output (planar_code), and counts separating pentagons over fullerene isomers. |
| Fullerene predicates & small test cages | `src/search_fullerene_localised.py`, `src/build_small_fullerenes.py` | `fullerene_ok`, `is_ipr`, planar `faces_of`, and a drum construction used for unit tests. |

### Reported totals (for cross-checking)

The acyclic-case run examines **6,623** pairs `(T, M)`, **917,415** maps `ψ`,
and **12,005,412** attachment patterns, yielding **12,183** complete
configuration records. Of these, **1,266** are eliminated by the exit
completion (A4) and **2,829** by the girth condition in (A5); each of the
remaining **8,088** is eliminated by the face conditions in (A5). **No
configuration record survives.**

## Prerequisites

- Python ≥ 3.9 and `pip install -r requirements.txt` (NetworkX).
- `gentreeg` (from the `nauty`/`Traces` suite) on `PATH` or in `./` — only for `gen_forests_nauty.py`.
- `buckygen` — only for `locate.py` / `sample.py`.

## Quick start: run everything

```bash
pip install -r requirements.txt
./run_all.sh          # ~2 min on one core; runs all self-contained checks
```
`run_all.sh` runs the three acyclic-case searches (k = 1, 2, 3), the
`gentreeg` cross-check (skipped with a note if nauty is absent), and the
explicit girth-5 witness verification, then prints a pass/fail summary.

## Reproducing the main finite verification

```bash
pip install -r requirements.txt
cd src
# Lemma (acyclic case): the three runs quoted in the paper.
python search_exit2.py 1 4 1 5     # k=1, c<=4, m in 1..5
python search_exit2.py 2 3 1 7     # k=2, c<=3, m in 1..7
python search_exit2.py 3 2 1 9     # k=3, c<=2, m in 1..9   (subsumes k=1,2)
```
Each run reports the per-filter elimination counts and ends with zero surviving
records. The whole verification runs in under two minutes on a single core.

## Reproducing the girth-5 witness

```bash
cd src
python verify_witness_paper.py   # verifies the EXACT 20-vertex graph from the paper
python witness.py                # heuristic search for a completion (optional)
```
`verify_witness_paper.py` confirms the paper's explicit graph is cubic,
3-connected, of girth 5, has trivial automorphism group, is non-planar, and
admits the stated rank-17 partial automorphism localised behind a 2-edge
interface. Expected output ends with `ALL CHECKS PASSED`.

## License

MIT — see [LICENSE](LICENSE).
