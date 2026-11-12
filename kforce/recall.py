"""Re-derive calls from stored kforce outputs with new thresholds, without re-analysis."""
from __future__ import annotations
import argparse
import csv
import glob
import os
from .config import decide_call

def recall_call(
    n_mut: int,
    vaf: float,
    n_mut_unique: int,
    evidence: str,
    min_mut_molecules: int,
    vaf_floor: float,
) -> str:
    """Return the call implied by stored counts under new thresholds."""
    return decide_call(n_mut, vaf, n_mut_unique, evidence, min_mut_molecules, vaf_floor)


def recall_row(row: dict, min_mut_molecules: int, vaf_floor: float) -> dict:
    """Update a result row's call and threshold columns in place."""
    row["call"] = recall_call(
        int(row["n_mut"]),
        float(row["vaf"]),
        int(row["n_mut_unique"]),
        row.get("evidence", "kmer"),
        min_mut_molecules,
        vaf_floor,
    )
    row["min_mut_molecules"] = min_mut_molecules
    row["vaf_floor"] = vaf_floor
    return row


def recall_tsv(in_path: str, out_path: str, mmm: int, floor: float) -> int:
    """Recall one results TSV and return the number of changed calls."""
    rows = list(csv.DictReader(open(in_path), delimiter="\t"))
    if not rows:
        return 0
    changed = 0
    cols = list(rows[0].keys())
    for r in rows:
        old = r["call"]
        recall_row(r, mmm, floor)
        changed += r["call"] != old
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return changed


def recall_dir(in_dir: str, out_dir: str, mmm: int, floor: float):
    """Recall every results TSV in a directory."""
    os.makedirs(out_dir, exist_ok=True)
    n_files = n_changed = 0
    for f in sorted(glob.glob(os.path.join(in_dir, "*.kforce.tsv"))):
        n_changed += recall_tsv(f, os.path.join(out_dir, os.path.basename(f)), mmm, floor)
        n_files += 1
    return n_files, n_changed


def main(argv=None) -> int:
    """Command-line entry point for recall."""
    ap = argparse.ArgumentParser(prog="kforce recall", description=__doc__)
    ap.add_argument(
        "--in", dest="inp", required=True, help="kforce.tsv file or a directory of them"
    )
    ap.add_argument("--out", required=True, help="output dir (recalled TSVs)")
    ap.add_argument("--min-mut-molecules", dest="mmm", type=int, required=True)
    ap.add_argument("--vaf-floor", dest="floor", type=float, default=0.0)
    args = ap.parse_args(argv)
    if os.path.isdir(args.inp):
        nf, nc = recall_dir(args.inp, args.out, args.mmm, args.floor)
        print(
            f"[recall] {nf} files, {nc} calls changed (mmm={args.mmm}, vaf_floor={args.floor}) "
            f"-> {args.out}"
        )
    else:
        os.makedirs(args.out, exist_ok=True)
        out = os.path.join(args.out, os.path.basename(args.inp))
        nc = recall_tsv(args.inp, out, args.mmm, args.floor)
        print(f"[recall] {nc} calls changed -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
