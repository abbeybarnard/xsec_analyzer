#!/usr/bin/env python3

import argparse
from pathlib import Path
import numpy as np
import uproot

# ---------------------------
# CONFIG
# ---------------------------
DEFAULT_FILES_LIST = Path("configs/files_to_process_analysis.txt")
TREE_PATH = "XSecAnalyzer/NuMICC1eNp"

BRANCHES = [
    "is_true_signal",
    "pandora_sel_nue_cc",
    "wc_sel_nue_cc",
    "sel_nue_cc_union",
]

# ---------------------------
# COUNTERS
# ---------------------------
class Counters:
    def __init__(self):
        self.total_true_signal = 0

        self.pandora_selected = 0
        self.pandora_selected_true_signal = 0

        self.wc_selected = 0
        self.wc_selected_true_signal = 0

        self.union_selected = 0
        self.union_selected_true_signal = 0

# ---------------------------
# UTIL
# ---------------------------
def read_file_list(path):
    files = []
    with open(path, "r") as f:
        for line in f:
            line = line.split("#")[0].strip()
            if not line:
                continue
            parts = line.split()
            files.append(Path(parts[0]))
    return files

def print_metrics(name, selected, selected_true, total_true):
    eff = selected_true / total_true if total_true > 0 else 0
    pur = selected_true / selected if selected > 0 else 0

    print(name)
    print(f"  selected             = {selected}")
    print(f"  selected true signal = {selected_true}")
    print(f"  efficiency           = {eff:.6f}")
    print(f"  purity               = {pur:.6f}")
    print()

# ---------------------------
# MAIN
# ---------------------------
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--files-list", default=str(DEFAULT_FILES_LIST))
    args = parser.parse_args()

    files = read_file_list(args.files_list)

    counters = Counters()

    for i, fpath in enumerate(files):
        print(f"[{i+1}/{len(files)}] {fpath}")

        try:
            with uproot.open(fpath) as f:
                tree = f[TREE_PATH]

                arrays = tree.arrays(BRANCHES, library="np")

                true_signal = arrays["is_true_signal"].astype(bool)
                pandora = arrays["pandora_sel_nue_cc"].astype(bool)
                wc = arrays["wc_sel_nue_cc"].astype(bool)
                union = arrays["sel_nue_cc_union"].astype(bool)

                counters.total_true_signal += true_signal.sum()

                counters.pandora_selected += pandora.sum()
                counters.pandora_selected_true_signal += (pandora & true_signal).sum()

                counters.wc_selected += wc.sum()
                counters.wc_selected_true_signal += (wc & true_signal).sum()

                counters.union_selected += union.sum()
                counters.union_selected_true_signal += (union & true_signal).sum()

        except Exception as e:
            print(f"  Skipping file: {e}")

    print("\n==============================")
    print("   FINAL SELECTION METRICS")
    print("==============================\n")

    print(f"Total true signal = {counters.total_true_signal}\n")

    print_metrics("Pandora selection",
                  counters.pandora_selected,
                  counters.pandora_selected_true_signal,
                  counters.total_true_signal)

    print_metrics("WireCell selection",
                  counters.wc_selected,
                  counters.wc_selected_true_signal,
                  counters.total_true_signal)

    print_metrics("Union selection",
                  counters.union_selected,
                  counters.union_selected_true_signal,
                  counters.total_true_signal)


if __name__ == "__main__":
    main()