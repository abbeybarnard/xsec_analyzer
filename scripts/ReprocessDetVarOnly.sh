#!/bin/bash
# Reprocess only the merge_det_* NuMI RHC DetVar files (the ones that use
# the "_cv"/"_det" merged-production tree convention) through the full
# compute_bdt_scores.py -> AddWeights.sh -> ReprocessNTuples.sh pipeline,
# without re-touching the already-correct nueMC/numuMC/extBNB/dirtMC/detVarCV
# outputs.
#
# AddWeights.sh and ReprocessNTuples.sh always truncate and rewrite their
# master list files (configs/files_to_process_bdt_weighted.txt and
# configs/files_to_process_analysis.txt) to contain only the entries they
# were just given. This wrapper backs those files up first and merges the
# freshly-written DetVar lines back into them afterward, since the output
# file paths are deterministic (based on the input file's basename) and stay
# valid across reruns.
#
# Usage: ./scripts/ReprocessDetVarOnly.sh
# Expects setup_xsec_analyzer.sh and the bdt-env virtualenv to already be
# sourced (as in run_all.sh).

set -e

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "${repo_dir}"

DETVAR_INPUT_LIST="configs/files_to_process_detvar_only.txt"
SCRATCH_DIR=$(mktemp -d)
BDT_LIST="${SCRATCH_DIR}/files_to_process_bdt_detvar_only.txt"

MASTER_BDT_WEIGHTED="configs/files_to_process_bdt_weighted.txt"
MASTER_ANALYSIS="configs/files_to_process_analysis.txt"
BACKUP_BDT_WEIGHTED="${SCRATCH_DIR}/files_to_process_bdt_weighted.backup.txt"
BACKUP_ANALYSIS="${SCRATCH_DIR}/files_to_process_analysis.backup.txt"

touch "${MASTER_BDT_WEIGHTED}" "${MASTER_ANALYSIS}"
cp "${MASTER_BDT_WEIGHTED}" "${BACKUP_BDT_WEIGHTED}"
cp "${MASTER_ANALYSIS}" "${BACKUP_ANALYSIS}"

echo "============================================"
echo "Stage 1: compute_bdt_scores.py (DetVar-only)"
echo "============================================"
python3 scripts/compute_bdt_scores.py \
  --files-list "${DETVAR_INPUT_LIST}" \
  --written-list "${BDT_LIST}"

echo ""
echo "============================================"
echo "Stage 2: AddWeights.sh (DetVar-only)"
echo "============================================"
./scripts/AddWeights.sh "${BDT_LIST}"

echo ""
echo "============================================"
echo "Stage 3: ReprocessNTuples.sh (DetVar-only)"
echo "============================================"
./scripts/ReprocessNTuples.sh \
  /exp/uboone/data/users/abarnard/analysis/selection_output \
  NuMICC1eNp \
  "${MASTER_BDT_WEIGHTED}"

echo ""
echo "============================================"
echo "Restoring full master file lists"
echo "============================================"

merge_lists() {
  local backup="$1" fresh="$2" out="$3"
  python3 - "${backup}" "${fresh}" "${out}" <<'PYEOF'
import sys
backup, fresh, out = sys.argv[1:4]
with open(fresh) as f:
    new_lines = [l for l in f if l.strip()]
new_types = {l.split()[1] for l in new_lines}
with open(backup) as f:
    kept = [l for l in f if l.strip() and l.split()[1] not in new_types]
with open(out, "w") as f:
    f.writelines(kept)
    f.writelines(new_lines)
PYEOF
}

merge_lists "${BACKUP_BDT_WEIGHTED}" "${MASTER_BDT_WEIGHTED}" "${SCRATCH_DIR}/merged_bdt_weighted.txt"
mv "${SCRATCH_DIR}/merged_bdt_weighted.txt" "${MASTER_BDT_WEIGHTED}"

merge_lists "${BACKUP_ANALYSIS}" "${MASTER_ANALYSIS}" "${SCRATCH_DIR}/merged_analysis.txt"
mv "${SCRATCH_DIR}/merged_analysis.txt" "${MASTER_ANALYSIS}"

echo "Done. Master lists restored with DetVar entries refreshed:"
echo "  ${MASTER_BDT_WEIGHTED}"
echo "  ${MASTER_ANALYSIS}"
echo "(scratch files kept at ${SCRATCH_DIR} for inspection)"
