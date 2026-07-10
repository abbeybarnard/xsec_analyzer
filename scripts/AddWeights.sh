#!/bin/bash
# =============================================================
# AddWeights.sh
#
# Pre-processing script to add beamline geometry weights and
# fake weights to NuMI BDT files before running ProcessNTuples.
#
# Reads:  configs/files_to_process_bdt.txt
# Writes: configs/files_to_process_bdt_weighted.txt
#         Weighted ROOT files to OUTPUT_DIR
#
# File type logic:
#   nueMC / numuMC  → AddBeamlineGeometryWeights
#   dirtMC          → AddBeamlineGeometryWeights, then AddFakeWeights
#                     (using first numuMC weighted file as reference)
#   extBNB          → copied as-is (data, no weights needed)
#   nueDV etc.      → AddBeamlineGeometryWeights (treated as MC overlay)
#
# Usage: ./scripts/AddWeights.sh [INPUT_FILE_LIST]
#   INPUT_FILE_LIST defaults to configs/files_to_process_bdt.txt
# =============================================================

set -e

# ── Configuration ─────────────────────────────────────────────
INPUT_LIST="${1:-configs/files_to_process_bdt.txt}"
OUTPUT_DIR="/exp/uboone/data/users/abarnard/analysis/bdt_files_weighted"
OUTPUT_LIST="configs/files_to_process_bdt_weighted.txt"

# ── Sanity checks ─────────────────────────────────────────────
if [ ! -f "${INPUT_LIST}" ]; then
    echo "ERROR: Input file list not found: ${INPUT_LIST}"
    exit 1
fi

if ! command -v AddBeamlineGeometryWeights &> /dev/null; then
    echo "ERROR: AddBeamlineGeometryWeights not found on PATH."
    echo "       Have you run setup_xsec_analyzer.sh?"
    exit 1
fi

if ! command -v AddFakeWeights &> /dev/null; then
    echo "ERROR: AddFakeWeights not found on PATH."
    echo "       Have you run setup_xsec_analyzer.sh?"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

# ── Helper: derive output filename ────────────────────────────
# Strips the input directory, prepends OUTPUT_DIR, adds _weighted
make_output_path() {
    local input_path="$1"
    local basename
    basename=$(basename "${input_path}" .root)
    echo "${OUTPUT_DIR}/${basename}_weighted.root"
}

# ── Pass 1: find first numuMC weighted file (needed as reference
#            for AddFakeWeights on dirtMC files) ───────────────
REFERENCE_WEIGHTS_FILE=""

echo "============================================"
echo "Pass 1: Finding numuMC reference file..."
echo "============================================"

while IFS= read -r line || [ -n "$line" ]; do
    [[ -z "$line" || "$line" == \#* ]] && continue
    read -r file_path file_type beam_mode extra_cols <<< "$line"
    if [[ "$file_type" == "numuMC" ]]; then
        REFERENCE_WEIGHTS_FILE=$(make_output_path "${file_path}")
        echo "  Reference weights file will be: ${REFERENCE_WEIGHTS_FILE}"
        break
    fi
done < "${INPUT_LIST}"

if [[ -z "${REFERENCE_WEIGHTS_FILE}" ]]; then
    echo "WARNING: No numuMC file found in ${INPUT_LIST}."
    echo "         dirtMC files cannot be processed without a reference."
    echo "         If you have no dirtMC files this is fine."
fi

# ── Pass 2: process all files ─────────────────────────────────
echo ""
echo "============================================"
echo "Pass 2: Processing all files..."
echo "============================================"

# Clear output list
> "${OUTPUT_LIST}"

LINE_NUM=0

while IFS= read -r line || [ -n "$line" ]; do
    # preserve blank lines and comments in output list
    if [[ -z "$line" || "$line" == \#* ]]; then
        echo "$line" >> "${OUTPUT_LIST}"
        continue
    fi

    LINE_NUM=$((LINE_NUM + 1))

    # parse all columns — preserve anything beyond the first three
    read -r file_path file_type beam_mode extra_cols <<< "$line"

    OUTPUT_FILE=$(make_output_path "${file_path}")

    # rebuild any extra columns for the output list
    if [[ -n "${extra_cols}" ]]; then
        extra_suffix=" ${extra_cols}"
    else
        extra_suffix=""
    fi

    echo ""
    echo "--------------------------------------------"
    echo "File ${LINE_NUM}: ${file_path}"
    echo "  Type:      ${file_type}"
    echo "  Beam mode: ${beam_mode}"
    echo "  Output:    ${OUTPUT_FILE}"
    echo "--------------------------------------------"

    case "${file_type}" in

        nueMC|numuMC|nueDV)
            echo "  → Running AddBeamlineGeometryWeights..."
            AddBeamlineGeometryWeights \
                "${file_path}" \
                "${beam_mode}" \
                "${OUTPUT_FILE}"
            echo "  ✓ Done"
            ;;

        dirtMC)
            if [[ -z "${REFERENCE_WEIGHTS_FILE}" ]]; then
                echo "  ERROR: Cannot process dirtMC without a numuMC reference file."
                exit 1
            fi
            if [[ ! -f "${REFERENCE_WEIGHTS_FILE}" ]]; then
                echo "  ERROR: Reference weights file not found: ${REFERENCE_WEIGHTS_FILE}"
                echo "         numuMC must be processed before dirtMC."
                exit 1
            fi
            # Step 1: add beamline geometry weights to dirt
            DIRT_INTERMEDIATE="${OUTPUT_DIR}/tmp_dirt_beamline_$(basename ${file_path})"
            echo "  → Step 1: Running AddBeamlineGeometryWeights on dirtMC..."
            AddBeamlineGeometryWeights \
                "${file_path}" \
                "${beam_mode}" \
                "${DIRT_INTERMEDIATE}"
            echo "  ✓ Beamline weights added"
            # Step 2: add fake systematic weights map using numuMC as reference
            echo "  → Step 2: Running AddFakeWeights (reference: ${REFERENCE_WEIGHTS_FILE})..."
            AddFakeWeights \
                "${DIRT_INTERMEDIATE}" \
                "${REFERENCE_WEIGHTS_FILE}" \
                "${OUTPUT_FILE}"
            # clean up intermediate file
            rm -f "${DIRT_INTERMEDIATE}"
            echo "  ✓ Done"
            ;;

        extBNB|data|ext*)
            echo "  → Data file — copying as-is (no weights needed)..."
            cp "${file_path}" "${OUTPUT_FILE}"
            echo "  ✓ Done"
            ;;

        detVar*)
            echo "  → detVar file — copying as-is (no weights needed; DV-type"
            echo "    systematics compare selected yields against detVarCV,"
            echo "    they don't use the weight_* branches added here)..."
            cp "${file_path}" "${OUTPUT_FILE}"
            echo "  ✓ Done"
            ;;

        *)
            echo "  WARNING: Unknown file type '${file_type}' — copying as-is."
            cp "${file_path}" "${OUTPUT_FILE}"
            ;;

    esac

    # append to output list, preserving all columns
    echo "${OUTPUT_FILE} ${file_type} ${beam_mode}${extra_suffix}" >> "${OUTPUT_LIST}"

done < "${INPUT_LIST}"

# ── Summary ───────────────────────────────────────────────────
echo ""
echo "============================================"
echo "All files processed successfully!"
echo ""
echo "Output files:     ${OUTPUT_DIR}/"
echo "New file list:    ${OUTPUT_LIST}"
echo "============================================"
cat "${OUTPUT_LIST}"