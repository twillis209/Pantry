#!/usr/bin/env python3
"""Assemble the ordered FUSION weight list from batch status tables.

If no phenotype produced weights (every one was skipped by FUSION, or the
modality has no phenotypes), FUSION.profile_wgt.R cannot run on the empty
list. With --profile and --profile-err, this script then writes those two
files itself: the profile's header line, and a note with the skip reasons.
"""

import argparse
import csv
from pathlib import Path

# Header of FUSION.profile_wgt.R's output (fusion_twas/utils/FUSION.profile_wgt.R).
PROFILE_MODELS = ["top1", "blup", "enet", "bslmm", "lasso"]
PROFILE_COLUMNS = (
    ["id", "nsnps", "hsq", "hsq.se", "hsq.pv"]
    + [f"{model}.r2" for model in PROFILE_MODELS]
    + [f"{model}.pv" for model in PROFILE_MODELS]
)


def parse_args():
    parser = argparse.ArgumentParser()
    # Zero status files when the modality has no phenotypes (no TWAS batches).
    parser.add_argument("--status", required=True, nargs="*", type=Path)
    parser.add_argument("--weights-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--profile", type=Path, help="Written only if there are no weights")
    parser.add_argument("--profile-err", type=Path, help="Written only if there are no weights")
    return parser.parse_args()


def read_status(status_paths):
    """(successful weight paths in BED order, skipped (phenotype_id, diagnostic) rows)."""
    rows = []
    skipped = []
    for status_path in status_paths:
        with status_path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["outcome"] == "success":
                    rows.append((int(row["bed_index"]), Path(row["weight_path"])))
                elif row["outcome"] == "skipped":
                    skipped.append((int(row["bed_index"]), row["phenotype_id"], row["diagnostic"]))
                else:
                    raise ValueError(f"Unexpected outcome in {status_path}: {row['outcome']}")
    rows.sort(key=lambda item: item[0])
    indices = [item[0] for item in rows]
    if len(indices) != len(set(indices)):
        raise ValueError("Duplicate BED indices in TWAS batch status tables")
    skipped.sort(key=lambda item: item[0])
    return [path for _, path in rows], [(pid, diag) for _, pid, diag in skipped]


def successful_weights(status_paths):
    return read_status(status_paths)[0]


def no_weights_note(n_status_files, skipped):
    if n_status_files == 0:
        return "No TWAS weights: the modality has no phenotypes"
    return f"No TWAS weights: FUSION skipped all {len(skipped)} phenotypes"


def write_empty_profile(profile, profile_err, note, skipped):
    """Stand-ins for FUSION.profile_wgt.R's output when there are no weights."""
    profile.write_text("\t".join(PROFILE_COLUMNS) + "\n")
    with profile_err.open("w") as handle:
        handle.write(f"{note} (FUSION.profile_wgt.R not run).\n")
        if skipped:
            handle.write("\nSkipped phenotypes (from the TWAS batch status tables):\n")
            handle.write("phenotype_id\tdiagnostic\n")
            for phenotype_id, diagnostic in skipped:
                handle.write(f"{phenotype_id}\t{diagnostic}\n")


def main():
    args = parse_args()
    paths, skipped = read_status(args.status)
    weights_dir = args.weights_dir.resolve()
    with args.output.open("w") as handle:
        for path in paths:
            handle.write(f"{path.resolve().relative_to(weights_dir)}\n")
    if not paths and args.profile is not None and args.profile_err is not None:
        note = no_weights_note(len(args.status), skipped)
        write_empty_profile(args.profile, args.profile_err, note, skipped)
        print(f"NOTE: {note}; wrote an empty weight list and profile")


if __name__ == "__main__":
    main()
