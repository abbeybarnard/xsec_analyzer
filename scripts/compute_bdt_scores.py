#!/usr/bin/env python3
"""
Preprocess ntuples to add XGBoost BDT scores using PyROOT.

Usage:

Single-file mode:
    python3 scripts/compute_bdt_scores.py input.root output.root --beam-mode FHC

Batch mode:
    python3 scripts/compute_bdt_scores.py

Batch mode reads `files_to_process.txt` and writes the created output file
paths to `files_to_process_bdt.txt`.

Expected `files_to_process.txt` formats per line:
    input.root file_type horn_current
    input.root file_type horn_current output.root

Blank lines and lines starting with # are ignored.

This script:
1. Loads FHC and RHC XGBoost BDT models (.model files)
2. Extracts the 13 input variables from the Pandora tree
3. Computes BDT scores for each event
4. Writes a new ROOT file containing:
   - The Pandora tree with a new 'bdt_score' branch
   - The three WireCell trees, unchanged
   All directory structure is preserved. TTree.CloneTree() is used throughout,
   so all branch types (including nested-jagged) are handled natively by ROOT.
"""

import sys
import array
import argparse
from pathlib import Path

import numpy as np
import xgboost as xgb
import ROOT

# Suppress ROOT info/warning banners; keep errors.
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError

REPO_DIR = Path(__file__).resolve().parent.parent

# Paths to the three WireCell trees to retain in the output.
DEFAULT_WIRECELL_TREES = [
    "wcpselection/T_eval",
    "wcpselection/T_KINEvars",
    "wcpselection/T_PFeval",
]

# Features that are simple per-event scalars (not vector-indexed).
_SCALAR_FEATURES = {
    "shr_tkfit_dedx_Y",
    "shr_tkfit_gap10_dedx_Y",
    "shr_tkfit_2cm_dedx_Y",
    "shrmoliereavg",
    "tksh_distance",
}

# Features indexed by shr_id (leading shower) rather than trk_id (leading track).
# Mirrors the notebook's compute_pfng2shrfrac / compute_pfng2shravrg functions
# which use shr_id - 1 as the index.
_SHR_INDEXED_FEATURES = {
    "pfng2shrfrac",
    "pfng2shravrg",
}

# ---------------------------------------------------------------------------
# Sentinel / bad-value cleaning
# Mirrors the notebook's replace_bad_with_nan() logic exactly:
#   1. +-inf        -> NaN
#   2. |x| > 0.1 * float32_max (~3.4e37) -> NaN   (float32 overflow junk)
#   3. |x| >= 9999  -> NaN                         (explicit reco-fail sentinels)
#   4. x == -1      -> NaN   (only for features in _SENTINEL_MINUS1_FEATURES,
#                             where -1 is physically impossible, e.g. distances)
# ---------------------------------------------------------------------------
_SENTINEL_THRESHOLD = 9999.0
_FLOAT32_OVERFLOW_THRESHOLD = 0.1 * np.finfo(np.float32).max  # ~3.4e37

# Features for which -1 is a reco-failure sentinel, not a valid value.
_SENTINEL_MINUS1_FEATURES = {"trkshrhitdist2"}


def _clean_features(features: np.ndarray, feature_order: list) -> np.ndarray:
    """
    Replace sentinel values, +-inf, and float32 overflows with NaN in-place,
    using the same rules as the notebook's replace_bad_with_nan().

    Parameters
    ----------
    features     : float32 array of shape (n_events, n_features)
    feature_order: list of feature name strings, one per column

    Returns
    -------
    The same array with bad values replaced by NaN.
    """
    for col_i, fname in enumerate(feature_order):
        col = features[:, col_i]

        # 1. +-inf -> NaN
        col = np.where(np.isinf(col), np.nan, col)

        # 2. float32 overflow-scale junk -> NaN
        col = np.where(np.abs(col) > _FLOAT32_OVERFLOW_THRESHOLD, np.nan, col)

        # 3. explicit reco-fail sentinels -> NaN
        col = np.where(np.abs(col) >= _SENTINEL_THRESHOLD, np.nan, col)

        # 4. -1 sentinel for specific features -> NaN
        if fname in _SENTINEL_MINUS1_FEATURES:
            col = np.where(col == -1.0, np.nan, col)

        features[:, col_i] = col

    return features


# ---------------------------------------------------------------------------
# Model loading
# ---------------------------------------------------------------------------

def load_models(configs_dir=None):
    """Load FHC and RHC XGBoost Booster models from disk."""
    if configs_dir is None:
        configs_dir = REPO_DIR / "configs" / "bdt"
    configs_dir = Path(configs_dir)
    if not configs_dir.is_absolute():
        configs_dir = REPO_DIR / configs_dir

    fhc_path = configs_dir / "fhc_bdt.model"
    rhc_path = configs_dir / "rhc_bdt.model"

    for p in (fhc_path, rhc_path):
        if not p.exists():
            raise FileNotFoundError(f"Model not found: {p}")

    print(f"Loading FHC model: {fhc_path}")
    fhc = xgb.Booster()
    fhc.load_model(str(fhc_path))

    print(f"Loading RHC model: {rhc_path}")
    rhc = xgb.Booster()
    rhc.load_model(str(rhc_path))

    return fhc, rhc


# ---------------------------------------------------------------------------
# Feature extraction
# ---------------------------------------------------------------------------

def _resolve_branch(tree, base_name):
    """
    Return the actual branch name in *tree*, trying plain name then
    the 'pandora_' prefix. Raises KeyError if neither exists.
    """
    if tree.GetBranch(base_name):
        return base_name
    prefixed = f"pandora_{base_name}"
    if tree.GetBranch(prefixed):
        return prefixed
    raise KeyError(
        f"Required branch not found: '{base_name}' (also tried '{prefixed}')"
    )


def _as_float32(value):
    """Safely cast a scalar or 0-d array-like to Python float."""
    try:
        return float(value)
    except TypeError:
        return float(np.asarray(value).flat[0])


def extract_bdt_features(tree, feature_order, max_events=None):
    """
    Extract BDT input features from a PyROOT TTree.

    Uses ROOT.RDataFrame.AsNumpy() for bulk extraction, then applies
    track-index and shower-index logic for vector branches.

    Indexing rules (matching the notebook exactly):
      - trk_id is one-based; trk_id=0 means no leading track -> NaN
      - shr_id is one-based; shr_id=0 means no leading shower -> NaN
      - pfng2hipfrac, pfng2hipavrg -> indexed by trk_id - 1
      - pfng2shrfrac, pfng2shravrg -> indexed by shr_id - 1
      - scalar features             -> read directly, no indexing
      - subcluster                  -> shrsubclusters0 + 1 + 2

    Returns
    -------
    np.ndarray of shape (n_events, n_features), dtype float32
    """
    aux_names = [
        "shrsubclusters0", "shrsubclusters1", "shrsubclusters2",
        "trk_id", "shr_id",
    ]

    all_base = []
    seen = set()
    for name in feature_order:
        if name != "subcluster" and name not in seen:
            all_base.append(name)
            seen.add(name)
    for name in aux_names:
        if name not in seen:
            all_base.append(name)
            seen.add(name)

    branch_map = {name: _resolve_branch(tree, name) for name in all_base}
    actual_branches = list(branch_map.values())

    rdf = ROOT.RDataFrame(tree)
    if max_events is not None:
        rdf = rdf.Range(int(max_events))

    print(f"  Reading {len(actual_branches)} branches with RDataFrame ...")
    np_dict = rdf.AsNumpy(columns=actual_branches)

    n_events = len(np_dict[branch_map["trk_id"]])
    if n_events == 0:
        raise ValueError("No events found while extracting features.")

    # trk_id: one-based, 0 = no leading track -> trk_idx = -1 (sentinel)
    raw_trk = np_dict[branch_map["trk_id"]]
    if hasattr(raw_trk[0], "__len__"):
        trk_idx = np.array(
            [int(v[0]) - 1 if len(v) > 0 else -1 for v in raw_trk],
            dtype=int,
        )
    else:
        trk_idx = raw_trk.astype(int) - 1
    # Do NOT clip -- leave -1 as sentinel for "no leading track".

    # shr_id: one-based, 0 = no leading shower -> shr_idx = -1 (sentinel)
    raw_shr = np_dict[branch_map["shr_id"]]
    if hasattr(raw_shr[0], "__len__"):
        shr_idx = np.array(
            [int(v[0]) - 1 if len(v) > 0 else -1 for v in raw_shr],
            dtype=int,
        )
    else:
        shr_idx = raw_shr.astype(int) - 1
    # Do NOT clip -- leave -1 as sentinel for "no leading shower".

    sub0 = np.asarray(np_dict[branch_map["shrsubclusters0"]], dtype=np.float32)
    sub1 = np.asarray(np_dict[branch_map["shrsubclusters1"]], dtype=np.float32)
    sub2 = np.asarray(np_dict[branch_map["shrsubclusters2"]], dtype=np.float32)
    subcluster = sub0 + sub1 + sub2

    features = np.zeros((n_events, len(feature_order)), dtype=np.float32)

    for col_i, fname in enumerate(feature_order):
        if fname == "subcluster":
            features[:, col_i] = subcluster
            continue

        actual = branch_map[fname]
        arr = np_dict[actual]

        if fname in _SCALAR_FEATURES:
            features[:, col_i] = np.asarray(arr, dtype=np.float32)
        else:
            col = np.zeros(n_events, dtype=np.float32)
            for j, v in enumerate(arr):
                # Choose the correct index for this feature:
                #   shower features -> shr_id - 1
                #   track features  -> trk_id - 1
                if fname in _SHR_INDEXED_FEATURES:
                    py_idx = int(shr_idx[j])
                else:
                    py_idx = int(trk_idx[j])

                # -1 means no leading track/shower -> NaN (missing value for XGBoost)
                if py_idx < 0:
                    col[j] = np.nan
                elif hasattr(v, "__len__"):
                    if len(v) > py_idx:
                        col[j] = _as_float32(v[py_idx])
                    elif len(v) > 0:
                        col[j] = _as_float32(v[0])
                else:
                    col[j] = _as_float32(v)
            features[:, col_i] = col

    return features


# ---------------------------------------------------------------------------
# Tree discovery
# ---------------------------------------------------------------------------

def find_pandora_tree(tfile):
    """
    Locate the Pandora / NeutrinoSelectionFilter tree in *tfile*.

    Tries the canonical path first. Some DetVar productions ("merge_det_..."
    files) instead store every tree twice, suffixed "_cv" (a paired
    central-value subset) and "_det" (the actual variation sample) -- for
    those we explicitly prefer "_det", since this analysis gets its detVar
    CV from the standalone detVarCV sample instead (see
    SystematicsCalculator.cxx, where the alternate-CV logic is gated on
    `!useNuMI`). Falls back to walking all TTrees and picking the one with
    the most required branches only if neither convention matches.

    Returns
    -------
    (path_str, TTree object, suffix_str)
    suffix_str is "_det" for merged-production files, "" otherwise. Sibling
    trees ending in this suffix should be kept (and renamed by stripping the
    suffix); their "_cv" counterparts should be dropped.
    """
    preferred = "nuselection/NeutrinoSelectionFilter"
    obj = tfile.Get(preferred)
    if obj and obj.InheritsFrom("TTree"):
        print(f"  Found Pandora tree at canonical path: {preferred}")
        return preferred, obj, ""

    det_path = preferred + "_det"
    obj = tfile.Get(det_path)
    if obj and obj.InheritsFrom("TTree"):
        print(f"  Found Pandora tree at merged-DetVar path: {det_path}")
        print(f"  (dropping paired '_cv' trees; using standalone detVarCV "
              f"sample for the CV instead)")
        return det_path, obj, "_det"

    required = {
        "shr_tkfit_dedx_Y",
        "shrmoliereavg",
        "tksh_distance",
        "trk_id",
        "shrsubclusters0",
    }

    best_path, best_tree, best_count = None, None, -1

    def walk(directory, prefix=""):
        nonlocal best_path, best_tree, best_count
        for key in directory.GetListOfKeys():
            name = key.GetName()
            classname = key.GetClassName()
            full = f"{prefix}/{name}" if prefix else name
            if classname.startswith("TTree"):
                t = key.ReadObj()
                branches = {b.GetName() for b in t.GetListOfBranches()}
                count = sum(
                    1 for r in required
                    if r in branches or f"pandora_{r}" in branches
                )
                if count > best_count:
                    best_count = count
                    best_path = full
                    best_tree = t
            elif classname.startswith("TDirectory"):
                walk(key.ReadObj(), full)

    walk(tfile)

    if best_path is None or best_count < len(required):
        raise ValueError(
            f"Could not find a TTree with all required Pandora branches. "
            f"Best: '{best_path}' ({best_count}/{len(required)} branches found)."
        )

    print(f"  Found Pandora tree (auto-detected): {best_path}")
    return best_path, best_tree, ""


# ---------------------------------------------------------------------------
# Output file helpers
# ---------------------------------------------------------------------------

def _get_or_mkdir(tfile, dir_path):
    """
    Navigate (or create) nested sub-directories inside *tfile*, returning
    the deepest TDirectory. An empty dir_path returns *tfile* itself.
    """
    current = tfile
    for part in dir_path.split("/"):
        if not part:
            continue
        existing = current.Get(part)
        if existing and existing.InheritsFrom("TDirectory"):
            current = existing
        else:
            new_dir = current.mkdir(part)
            if not new_dir:
                raise RuntimeError(
                    f"Failed to mkdir '{part}' inside {current.GetPath()}"
                )
            current = new_dir
    return current


def copy_tree_with_score(src_dir, tree_path, dst_dir, scores=None, max_events=None,
                          output_name=None):
    """
    Clone a TTree from src_dir into dst_dir, optionally renaming it to
    *output_name* (used to strip a merged-production "_det" suffix so the
    output file uses the canonical tree name expected downstream).

    dst_dir is already the correct directory -- do NOT recreate paths here.
    """
    tree_name = tree_path.split("/")[-1]
    if output_name is None:
        output_name = tree_name

    src_tree = src_dir.Get(tree_name)
    if not src_tree or not src_tree.InheritsFrom("TTree"):
        raise ValueError(f"TTree not found: '{tree_path}'")

    dst_dir.cd()

    # Case 1: no BDT score -> direct clone
    if scores is None:
        n_clone = -1 if max_events is None else int(max_events)
        clone = src_tree.CloneTree(n_clone)
        if clone is None:
            raise RuntimeError(f"CloneTree failed for '{tree_path}'")
        if output_name != tree_name:
            clone.SetName(output_name)
        clone.Write("", ROOT.TObject.kOverwrite)
        print(f"  Copied: {tree_path} -> {output_name} ({clone.GetEntries()} entries)")
        return clone

    # Case 2: add BDT score
    n_entries = src_tree.GetEntries() if max_events is None else min(
        int(max_events),
        int(src_tree.GetEntries()),
    )

    clone = src_tree.CloneTree(0)
    if clone is None:
        raise RuntimeError(f"CloneTree failed for '{tree_path}'")
    if output_name != tree_name:
        clone.SetName(output_name)

    bdt_val = array.array("f", [0.0])
    clone.Branch("bdt_score", bdt_val, "bdt_score/F")

    for i in range(n_entries):
        src_tree.GetEntry(i)
        bdt_val[0] = float(scores[i]) if i < len(scores) else 0.0
        clone.Fill()

        if i > 0 and i % 10000 == 0:
            print(f"  Filled {i} entries in {tree_path}")

    clone.Write("", ROOT.TObject.kOverwrite)
    print(f"  Copied: {tree_path} -> {output_name} (+ bdt_score)")
    return clone


def _copy_object_recursive(src_dir, dst_dir, target_tree_path, score, max_events,
                            suffix="", current_path=""):
    """
    Recursively copy only the Pandora and WireCell parts of the file.

    Rules:
      - keep only top-level directories: nuselection, wcpselection
      - inside those, copy all TTrees except wcpselection/T_spacepoints
      - add bdt_score only to the Pandora tree
      - drop everything else
      - if *suffix* is non-empty (merged CV/DetVar production files, see
        find_pandora_tree), only trees ending in *suffix* are kept, renamed
        by stripping the suffix; their paired "_cv" siblings are dropped.
    """
    for key in src_dir.GetListOfKeys():
        obj = key.ReadObj()
        name = obj.GetName()
        class_name = obj.ClassName() if hasattr(obj, 'ClassName') else type(obj).__name__
        full_path = f"{current_path}/{name}" if current_path else name

        print(f"[DEBUG] Found object: {full_path} (type: {class_name})")

        # Keep only desired top-level directories
        if current_path == "" and obj.InheritsFrom("TDirectory"):
            if name not in {"nuselection", "wcpselection"}:
                print(f"[DEBUG] Skipping top-level directory: {name}")
                continue

        if obj.InheritsFrom("TDirectory"):
            existing = dst_dir.Get(name)
            if existing and existing.InheritsFrom("TDirectory"):
                out_subdir = existing
            else:
                out_subdir = dst_dir.mkdir(name)
                if not out_subdir:
                    raise RuntimeError(
                        f"Failed to mkdir '{name}' inside {dst_dir.GetPath()}"
                    )

            _copy_object_recursive(
                obj,
                out_subdir,
                target_tree_path=target_tree_path,
                score=score,
                max_events=max_events,
                suffix=suffix,
                current_path=full_path,
            )

        elif obj.InheritsFrom("TTree"):
            out_name = name
            if suffix:
                if name.endswith("_cv"):
                    print(f"[DEBUG] Dropping paired-CV tree (standalone "
                          f"detVarCV sample supplies the CV instead): {full_path}")
                    continue
                if name.endswith(suffix):
                    out_name = name[: -len(suffix)]
                else:
                    print(f"[DEBUG] Warning: '{full_path}' has neither "
                          f"'{suffix}' nor '_cv' suffix; copying unchanged")

            # Skip specific tree (compare against the *renamed* output path)
            out_full_path = f"{current_path}/{out_name}" if current_path else out_name
            if out_full_path == "wcpselection/T_spacepoints":
                print(f"[DEBUG] Skipping tree: {full_path}")
                continue

            try:
                if full_path == target_tree_path:
                    copy_tree_with_score(
                        src_dir,
                        full_path,
                        dst_dir,
                        scores=score,
                        max_events=max_events,
                        output_name=out_name,
                    )
                else:
                    copy_tree_with_score(
                        src_dir,
                        full_path,
                        dst_dir,
                        scores=None,
                        max_events=max_events,
                        output_name=out_name,
                    )
            except ValueError as e:
                print(f"[DEBUG] Warning: {e}")

        else:
            print(f"[DEBUG] Skipping non-tree object: {full_path}")


# ---------------------------------------------------------------------------
# File list parsing and batch mode
# ---------------------------------------------------------------------------

_DEFAULT_BDT_OUTPUT_DIR = Path("/exp/uboone/data/users/abarnard/analysis/bdt_files")


def _resolve_path_maybe_relative(path_str, base_dir):
    p = Path(path_str)
    if p.is_absolute():
        return p
    return (base_dir / p).resolve()


def _default_bdt_output_file(input_file):
    """
    Given an input file path, return the default output file path in the BDT output directory.
    """
    input_file = Path(input_file)
    out_name = f"{input_file.stem}_bdt.root"
    return _DEFAULT_BDT_OUTPUT_DIR / out_name


def read_files_to_process(list_path):
    """
    Read files_to_process.txt.

    Supported line formats:
        input.root file_type horn_current
        input.root file_type horn_current output.root

    Returns a list of dicts with keys:
        input_file, file_type, beam_mode, output_file
    """
    list_path = Path(list_path)
    if not list_path.exists():
        raise FileNotFoundError(f"File list not found: {list_path}")

    base_dir = list_path.parent
    entries = []

    with list_path.open("r", encoding="utf-8") as f:
        for lineno, raw_line in enumerate(f, start=1):
            line = raw_line.split("#", 1)[0].strip()
            if not line:
                continue

            parts = line.replace(",", " ").split()
            if len(parts) < 3:
                raise ValueError(
                    f"{list_path}:{lineno}: expected at least "
                    f"'input.root file_type horn_current'"
                )

            input_file = _resolve_path_maybe_relative(parts[0], base_dir)
            file_type = parts[1]
            beam_mode = parts[2].upper()

            if beam_mode not in {"FHC", "RHC"}:
                raise ValueError(
                    f"{list_path}:{lineno}: horn current must be FHC or RHC, got '{parts[2]}'"
                )

            if len(parts) >= 4:
                output_file = _resolve_path_maybe_relative(parts[3], base_dir)
            else:
                output_file = _default_bdt_output_file(input_file)

            entries.append(
                {
                    "input_file": input_file,
                    "file_type": file_type,
                    "beam_mode": beam_mode,
                    "output_file": output_file,
                }
            )

    if not entries:
        raise ValueError(f"No usable entries found in {list_path}")

    return entries


def _output_already_processed(output_file):
    """
    Return True if output_file already contains a complete, valid BDT-scored
    Pandora tree (i.e. a previous run finished it successfully).

    Used to make batch mode resumable: if one entry in files_to_process.txt
    has a bad path and kills the run, re-running should skip everything that
    already succeeded rather than reprocessing it from scratch.
    """
    output_file = Path(output_file)
    if not output_file.exists() or output_file.stat().st_size == 0:
        return False

    tfile = ROOT.TFile.Open(str(output_file), "READ")
    try:
        if not tfile or tfile.IsZombie():
            return False
        try:
            _, tree, _ = find_pandora_tree(tfile)
        except (ValueError, KeyError):
            return False
        return bool(tree.GetBranch("bdt_score")) and tree.GetEntries() > 0
    finally:
        try:
            tfile.Close()
        except Exception:
            pass


def process_batch(
    files_list="files_to_process.txt",
    written_list="files_to_process_bdt.txt",
    max_events=None,
    config_dir=None,
):
    """
    Process every file listed in files_to_process.txt and write the output
    filenames to files_to_process_bdt.txt as:

        output.root file_type horn_current

    Entries whose output file already exists and is complete (see
    _output_already_processed) are skipped, so a failed/interrupted batch
    run can simply be re-run to pick up where it left off.
    """
    entries = read_files_to_process(files_list)
    written_list = Path(written_list)
    written_list.parent.mkdir(parents=True, exist_ok=True)

    with written_list.open("w", encoding="utf-8") as out:
        for idx, entry in enumerate(entries, start=1):
            input_file = entry["input_file"]
            file_type = entry["file_type"]
            beam_mode = entry["beam_mode"]
            output_file = entry["output_file"]

            if _output_already_processed(output_file):
                print(
                    f"\n[{idx}/{len(entries)}] {input_file} "
                    f"({file_type}, {beam_mode}) -> {output_file}"
                    " [already processed, skipping]"
                )
            else:
                print(
                    f"\n[{idx}/{len(entries)}] {input_file} "
                    f"({file_type}, {beam_mode}) -> {output_file}"
                )

                process_ntuple(
                    input_file,
                    output_file,
                    beam_mode=beam_mode,
                    max_events=max_events,
                    config_dir=config_dir,
                )

            out.write(f"{output_file} {file_type} {beam_mode}\n")
            out.flush()

    print(f"\nWrote output file list to: {written_list}")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def process_ntuple(
    input_file,
    output_file,
    beam_mode="FHC",
    max_events=None,
    config_dir=None,
):
    """
    Full pipeline:
      1. Load models.
      2. Open input file, locate Pandora tree.
      3. Extract features, compute BDT scores.
      4. Write output file with Pandora tree (+ bdt_score) and WireCell trees.
    """
    input_file = Path(input_file)
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    fhc_booster, rhc_booster = load_models(config_dir)
    booster = rhc_booster if beam_mode == "RHC" else fhc_booster

    feature_names = booster.feature_names or [
        "shr_tkfit_dedx_Y",
        "shr_tkfit_gap10_dedx_Y",
        "shr_tkfit_2cm_dedx_Y",
        "pfng2hipfrac",
        "pfng2hipavrg",
        "ng2hip_r1cm",
        "ng2hip_r10cm",
        "pfng2shrfrac",
        "pfng2shravrg",
        "shrmoliereavg",
        "subcluster",
        "tksh_distance",
        "trkshrhitdist2",
    ]

    print(f"\nInput  : {input_file}")
    print(f"Output : {output_file}")
    print(f"Mode   : {beam_mode}")

    # Open with PyROOT and locate the Pandora tree.
    src = ROOT.TFile.Open(str(input_file), "READ")
    if not src or src.IsZombie():
        raise IOError(f"Cannot open input file: {input_file}")

    missing = []
    try:
        print("\n[1/3] Locating Pandora tree ...")
        pandora_path, pandora_tree, pandora_suffix = find_pandora_tree(src)
        print(f"Using Pandora tree: {pandora_path}")

        print(f"\n[2/3] Extracting features and computing BDT scores ...")
        features = extract_bdt_features(
            pandora_tree,
            feature_order=feature_names,
            max_events=max_events,
        )
        n_events = features.shape[0]
        print(f"  Events : {n_events}, Features : {features.shape[1]}")

        print(f"  Cleaning sentinel / bad values (matching notebook logic) ...")
        features = _clean_features(features, feature_names)

        dmatrix = xgb.DMatrix(features, feature_names=feature_names,
                               missing=np.nan)
        scores = booster.predict(dmatrix)
        print(f"  Score range : [{scores.min():.4f}, {scores.max():.4f}]")

        print(f"\n[3/3] Writing output file ...")
        dst = ROOT.TFile.Open(str(output_file), "RECREATE")
        if not dst or dst.IsZombie():
            raise IOError(f"Cannot create output file: {output_file}")

        try:
            # Copy everything, adding bdt_score only to the Pandora tree.
            _copy_object_recursive(
                src,
                dst,
                target_tree_path=pandora_path,
                score=scores,
                max_events=max_events,
                suffix=pandora_suffix,
            )

            # Track missing WireCell trees for reporting.
            for wc_path in DEFAULT_WIRECELL_TREES:
                check_path = wc_path + pandora_suffix if pandora_suffix else wc_path
                obj = src.Get(check_path)
                if not (obj and obj.InheritsFrom("TTree")):
                    missing.append(wc_path)

            if missing:
                print("\n  Warning - WireCell trees not found in source (skipped):")
                for p in missing:
                    print(f"    {p}")

            dst.Write("", ROOT.TObject.kOverwrite)

        finally:
            try:
                dst.Close()
            except Exception:
                pass

    finally:
        try:
            src.Close()
        except Exception:
            pass

    print(f"\nDone. Output written to: {output_file}")
    print(f"  Pandora tree : {pandora_path}  (+ bdt_score)")
    kept = [p for p in DEFAULT_WIRECELL_TREES if p not in missing]
    print(f"  WireCell trees kept : {kept}")

    return str(output_file)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Add XGBoost BDT scores to ntuples using PyROOT",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input", nargs="?", help="Input ROOT ntuple file")
    parser.add_argument("output", nargs="?", help="Output ROOT ntuple file")
    parser.add_argument(
        "--beam-mode",
        default="FHC",
        choices=["FHC", "RHC"],
        help="Beam mode -- used only for single-file mode",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=None,
        help="Maximum number of events to process (default: all)",
    )
    parser.add_argument(
        "--config-dir",
        default=str(REPO_DIR / "configs" / "bdt"),
        help="Directory containing fhc_bdt.model / rhc_bdt.model",
    )
    parser.add_argument(
        "--files-list",
        default=str(REPO_DIR / "configs" / "files_to_process.txt"),
        help="Batch input list file (used when no positional input/output are given)",
    )
    parser.add_argument(
        "--written-list",
        default=str(REPO_DIR / "configs" / "files_to_process_bdt.txt"),
        help="Batch output list file containing the created ROOT files",
    )

    args = parser.parse_args()

    try:
        if (args.input is None) != (args.output is None):
            parser.error("Either provide both input and output, or provide neither.")

        if args.input is not None:
            output_file = args.output
            if output_file is None:
                output_file = str(_default_bdt_output_file(args.input))
            process_ntuple(
                args.input,
                output_file,
                beam_mode=args.beam_mode,
                max_events=args.max_events,
                config_dir=args.config_dir,
            )
        else:
            process_batch(
                files_list=args.files_list,
                written_list=args.written_list,
                max_events=args.max_events,
                config_dir=args.config_dir,
            )

    except Exception as exc:
        print(f"\nError: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()