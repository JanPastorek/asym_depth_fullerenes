#!/usr/bin/env bash
# Reproduce every self-contained finite verification in the paper
#   "Extremal Asymmetric Depth on Surfaces and Hidden Almost-Global
#    Symmetries in IPR Fullerenes".
#
# Usage:   ./run_all.sh
# Runtime: about two minutes on a single core.
# Requires: Python 3.9+ and networkx (pip install -r requirements.txt).
# Optional: `gentreeg` (nauty suite) on PATH enables the independent
#           trusted-generator cross-check; skipped with a note if absent.
set -u
cd "$(dirname "$0")/src"

PY="${PYTHON:-python3}"
fail=0
step() { printf '\n=====  %s  =====\n' "$1"; }

step "Environment"
"$PY" --version || { echo "python3 not found"; exit 2; }
"$PY" - <<'PYCHK' || { echo "ERROR: networkx missing -> pip install -r requirements.txt"; exit 2; }
import networkx; print("networkx", networkx.__version__)
PYCHK

# ---------------------------------------------------------------------------
step "Lemma (acyclic case): exhaustive configuration search, k in {1,2,3}"
# Paper's three runs: (k, cmax) in {(1,4),(2,3),(3,2)}, m in 1..3+2k.
run_exit() {  # args: k cmax mhi
  local k=$1 cmax=$2 mhi=$3
  echo "--- search_exit2.py k=$k c<=$cmax m=1..$mhi ---"
  if "$PY" search_exit2.py "$k" "$cmax" 1 "$mhi" | tee /tmp/_exit_$k.out | tail -n 8; then
    if grep -q "no survivors" /tmp/_exit_$k.out; then
      echo "  RESULT: zero survivors (as claimed)."
    else
      echo "  RESULT: !!! SURVIVORS FOUND -- investigate !!!"; fail=1
    fi
  else
    echo "  RESULT: run failed"; fail=1
  fi
}
run_exit 1 4 5
run_exit 2 3 7
run_exit 3 2 9

# ---------------------------------------------------------------------------
step "Independent cross-check: subcubic carrier trees via nauty gentreeg"
if command -v gentreeg >/dev/null 2>&1 || [ -x ./gentreeg ]; then
  "$PY" gen_forests_nauty.py && echo "  gentreeg tree counts printed above."
else
  echo "  SKIPPED: 'gentreeg' (nauty suite) not found on PATH or in ./ ."
  echo "  Install nauty and re-run to reproduce the trusted-generator check."
fi

# ---------------------------------------------------------------------------
step "Proposition (girth-5 witness): verify the explicit 20-vertex graph"
if "$PY" verify_witness_paper.py; then
  echo "  RESULT: witness verified."
else
  echo "  RESULT: !!! witness verification FAILED !!!"; fail=1
fi

# ---------------------------------------------------------------------------
step "Summary"
if [ "$fail" -eq 0 ]; then
  echo "ALL SELF-CONTAINED VERIFICATIONS PASSED."
  echo "(The buckygen separating-pentagon census in locate.py/sample.py needs"
  echo " the buckygen binary and is run separately; see README.)"
  exit 0
else
  echo "ONE OR MORE CHECKS FAILED -- see output above."
  exit 1
fi
