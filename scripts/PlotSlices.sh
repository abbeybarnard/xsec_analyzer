# Number of expected command-line arguments
num_expected=5

if [ "$#" -ne "$num_expected" ]; then
  echo "Usage: ./PlotSlices.sh FPM_Config SYST_Config SLICE_Config Univ_File PlotOutputDir"
  exit 1
fi

FPM_Config=$1
SYST_Config=$2
SLICE_Config=$3
Univ_File=$4
PlotOutputDir=$5

if [ ! -f "${FPM_Config}" ]; then
  echo "FPM_Config \"${FPM_Config}\" not found"
  exit 1
fi

if [ ! -f "${SYST_Config}" ]; then
  echo "SYST_Config \"${SYST_Config}\" not found"
  exit 2
fi

if [ ! -f "${SLICE_Config}" ]; then
  echo "SLICE_Config \"${SLICE_Config}\" not found"
  exit 3
fi

if [ ! -f "${Univ_File}" ]; then
  echo "Universe File \"${Univ_File}\" not found"
  exit 4
fi

if [ ! -d "${PlotOutputDir}" ]; then
  echo "Output directory \"${PlotOutputDir}\" not found"
  exit 5
fi

# SystematicsCalculator caches its POT-summed universes in a 'total_*'
# TDirectoryFile inside Univ_File the first time it's used, then silently
# reuses that cache forever after -- even once file_properties.txt, the
# systematics config, or the selection code have changed. Force a rebuild by
# default so this never serves stale numbers without you noticing. Set
# XSEC_ALLOW_CACHE=1 to skip this (e.g. for a final, already-validated
# large-scale run where you want the speed of reusing the existing cache).
if [ -z "${XSEC_ALLOW_CACHE:-}" ]; then
  SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  root -l -b -q "${SCRIPT_DIR}/invalidate_total_cache.C(\"${Univ_File}\")"
fi

SlicePlots ${FPM_Config} ${SYST_Config} ${SLICE_Config} ${Univ_File} ${PlotOutputDir}
