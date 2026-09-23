#!/usr/bin/env bash
# Check that the current solver reproduces an earlier revision bit-for-bit.
#
#     tools/verify_bitwise.sh [REF] [IDX] [TFIN] [TPRIN]
#
# Defaults: REF=HEAD~1, IDX=1, TFIN=20, TPRIN=10 — a couple of snapshots, not a
# full run. Both revisions get the same t_fin/t_prin override, so the only
# difference is the solver code itself.
#
# Run this on a GPU node. It needs enough device memory for the grid in DF.csv
# row IDX; for the published L=50 row that means a 40 GB-class card.
set -euo pipefail

REF=${1:-HEAD~1}
IDX=${2:-1}
TFIN=${3:-20}
TPRIN=${4:-10}

ROOT=$(git rev-parse --show-toplevel)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

echo "reference : $REF ($(git -C "$ROOT" rev-parse --short "$REF"))"
echo "current   : $(git -C "$ROOT" rev-parse --short HEAD)"
echo "sim index : $IDX   t_fin=$TFIN t_prin=$TPRIN"
echo

# Export both revisions.
mkdir -p "$WORK/ref" "$WORK/cur"
git -C "$ROOT" archive "$REF" | tar -x -C "$WORK/ref"
git -C "$ROOT" archive HEAD   | tar -x -C "$WORK/cur"

# Locate the solver directory in each tree. It is not hardcoded because the
# directory was renamed from FFT_2D_P_L50 to 2D, so older revisions differ.
solver_dir() {
    local found
    found=$(find "$1" -maxdepth 2 -name 2D.jl -print -quit)
    [ -n "$found" ] || { echo "no 2D.jl under $1" >&2; exit 1; }
    dirname "$found"
}

# Shorten the run identically in both copies.
for d in ref cur; do
    f="$(solver_dir "$WORK/$d")/InputParameters.jl"
    sed -i.bak -E "s/^t_fin = .*/t_fin = $TFIN/; s/^t_prin = .*/t_prin = $TPRIN/" "$f"
    rm -f "$f.bak"
done

for d in ref cur; do
    echo "--- running $d ---"
    ( cd "$(solver_dir "$WORK/$d")" \
      && DATA_DIR="$WORK/out_$d/" SIM_IDX="$IDX" \
         julia --project="$WORK/$d" --optimize=3 2D.jl )
done

echo
echo "--- comparing ---"
julia --project="$ROOT" "$ROOT/tools/compare_snapshots.jl" \
      "$WORK/out_ref/$IDX/Data" "$WORK/out_cur/$IDX/Data"
