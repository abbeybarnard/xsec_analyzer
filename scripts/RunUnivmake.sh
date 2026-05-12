#!/bin/bash
# =============================================================
# RunUnivmake.sh
#
# Runs univmake for each observable bin configuration.
# Each run takes ~1.5 hours, so run this overnight in a
# screen or tmux session:
#
#   screen -S univmake
#   source setup_xsec_analyzer.sh
#   ./scripts/RunUnivmake.sh v1_dummy
#
# Output ROOT files are written to OUTPUT_DIR/VERSION/.
#
# univmake argument order:
#   univmake LIST_FILE SYST_CONFIG OUTPUT_FILE FILE_PROPERTIES
#
# The bin config is NOT passed to univmake directly — it is
# read internally via the selection name in selections.conf.
# The LIST_FILE is files_to_process_analysis.txt.
# =============================================================

set -e

# ── Version argument ──────────────────────────────────────────
if [ "$#" -ne 1 ]; then
    echo "Usage: ./scripts/RunUnivmake.sh VERSION"
    echo "  e.g. ./scripts/RunUnivmake.sh v1_dummy"
    exit 1
fi

VERSION="$1"

# ── Configuration ─────────────────────────────────────────────
FILES="configs/files_to_process_analysis.txt"
SYST="configs/systcalc_numi.conf"
FILE_PROPS="${XSEC_ANALYZER_DIR}/configs/file_properties.txt"
OUTPUT_DIR="/exp/uboone/data/users/abarnard/analysis/univmake_output/${VERSION}"

echo "============================================"
echo "RunUnivmake version: ${VERSION}"
echo "Output directory:    ${OUTPUT_DIR}"
echo "============================================"

# ── Sanity checks ─────────────────────────────────────────────
if [ ! -f "${FILES}" ]; then
    echo "ERROR: File list not found: ${FILES}"
    exit 1
fi

if [ ! -f "${SYST}" ]; then
    echo "ERROR: Systematics config not found: ${SYST}"
    exit 1
fi

if [ ! -f "${FILE_PROPS}" ]; then
    echo "ERROR: File properties not found: ${FILE_PROPS}"
    exit 1
fi

if ! command -v univmake &> /dev/null; then
    echo "ERROR: univmake not found on PATH."
    echo "       Have you run setup_xsec_analyzer.sh?"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

# ── Helper: run univmake for one observable ───────────────────
run_univmake() {
    local label="$1"
    local bin_config="$2"
    local output_file="${OUTPUT_DIR}/univmake_${label}.root"

    if [ ! -f "${bin_config}" ]; then
        echo "ERROR: Bin config not found: ${bin_config}"
        exit 1
    fi

    echo ""
    echo "============================================"
    echo "Running univmake: ${label}"
    echo "  File list:       ${FILES}"
    echo "  Bin config:      ${bin_config}"
    echo "  Syst config:     ${SYST}"
    echo "  File properties: ${FILE_PROPS}"
    echo "  Output:          ${output_file}"
    echo "  Started:         $(date)"
    echo "============================================"

    univmake \
        "${FILES}" \
        "${bin_config}" \
        "${output_file}" \
        "${FILE_PROPS}"

    echo "  Finished:   $(date)"
    echo "  ✓ Done: ${output_file}"
}

# ── Run for each observable ───────────────────────────────────
# run_univmake "cos_opening_angle"   "configs/nuecc_cos_opening_angle_bin_config.txt"
# run_univmake "electron_energy"     "configs/nuecc_electron_energy_bin_config.txt"
# run_univmake "proton_ke"           "configs/nuecc_proton_ke_bin_config.txt"
run_univmake "double_differential" "configs/nuecc_double_differential_bin_config.txt"

# ── Summary ───────────────────────────────────────────────────
echo ""
echo "============================================"
echo "All univmake runs complete! [${VERSION}]"
echo "Output files:"
for f in "${OUTPUT_DIR}"/univmake_*.root; do
    echo "  ${f}"
done
echo "============================================"