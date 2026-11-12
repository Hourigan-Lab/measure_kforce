"""Command-line interface: kforce call | inspect | recall | aggregate."""
from __future__ import annotations
import argparse
import csv
import dataclasses
import glob
import sys
from pathlib import Path
import pysam
from .baseline import ConsensusSource, build_patient_model, write_baseline_detail
from .config import PRESETS, Params
from .inspect import write_inspection
from .manifest import filter_by_indel_size, load_manifest
from .pipeline import force_call, force_call_combined
from .reads import open_alignment
from .report import TIDY_COLUMNS, write_reads, write_report
from .support_bam import write_support_bam


def _build_params(args) -> Params:
    """Build run parameters from a preset plus command-line overrides."""
    base = PRESETS.get(args.level, PRESETS["wgs"])
    overrides = {}
    for field in (
        "k",
        "flank",
        "window",
        "max_mismatch",
        "min_mut_kmers",
        "min_wt_kmers",
        "min_mut_molecules",
        "vaf_floor",
        "max_reads_per_locus",
        "min_mapq",
        "baseline_haplotype_margin",
    ):
        val = getattr(args, field, None)
        if val is not None:
            overrides[field] = val
    if getattr(args, "bothend_parallel", False):
        overrides["bothend_parallel"] = True
    overrides["input_level"] = args.level
    return dataclasses.replace(base, **overrides)


def cmd_call(args) -> int:
    """Run the 'call' subcommand: force-call manifest variants in a BAM/CRAM."""
    params = _build_params(args)
    Path(args.out).mkdir(parents=True, exist_ok=True)
    sample = args.sample or Path(args.bam).name.split(".")[0]

    variants = load_manifest(args.manifest, patient=args.patient)
    if not variants:
        print("no variants loaded from manifest (check --patient filter)", file=sys.stderr)
        return 2

    variants, n_small = filter_by_indel_size(variants, args.min_indel_size)
    if n_small:
        print(
            f"[kforce] skipped {n_small} small indels (size <= {args.min_indel_size} bp) "
            f"per --min-indel-size"
        )
    if not variants:
        print("no variants remain after --min-indel-size filtering", file=sys.stderr)
        return 2

    if args.combined and not args.baseline_bam:
        print("--combined requires --baseline-bam", file=sys.stderr)
        return 2

    th = getattr(args, "threads", 1)
    fa = pysam.FastaFile(args.reference)
    bam = open_alignment(args.bam, reference=args.reference, threads=th)
    baseline_ref = getattr(args, "baseline_reference", None) or args.reference
    baseline_fa = pysam.FastaFile(baseline_ref) if baseline_ref != args.reference else fa
    baseline_bam = None
    details: list = []

    if args.combined:
        reference = "combined"
        baseline_bam = open_alignment(args.baseline_bam, reference=baseline_ref, threads=th)
        calls, refs, skipped, read_rows, details, support = force_call_combined(
            variants,
            fa,
            bam,
            params,
            baseline_bam,
            emit_reads=args.emit_reads,
            emit_bam=args.emit_bam,
            baseline_fa=baseline_fa,
        )
        rows = [c.as_row(params, sample, reference=r) for c, r in zip(calls, refs)]
    else:
        reference = "genome"
        model_builder = seq_source = None
        if args.baseline_bam:
            reference = "baseline"
            baseline_bam = open_alignment(args.baseline_bam, reference=baseline_ref, threads=th)
            seq_source = ConsensusSource(baseline_bam, baseline_fa, params)

            def model_builder(v):
                """Build the patient-specific model for one variant."""
                return build_patient_model(v, baseline_bam, baseline_fa, params)

        calls, skipped, read_rows, details, support = force_call(
            variants,
            fa,
            bam,
            params,
            emit_reads=args.emit_reads,
            model_builder=model_builder,
            seq_source=seq_source,
            emit_bam=args.emit_bam,
        )
        rows = [c.as_row(params, sample, reference=reference) for c in calls]
    if args.emit_bam:
        bam_path = str(Path(args.out) / f"{sample}.kforce.support.bam")
        info = write_support_bam(bam, support, bam_path, params)
        print(
            f"[kforce] support BAM: {info['n_mut_reads']} MUT + {info['n_wt_reads']} WT "
            f"reads (KF tag) over {info['n_variants']} variants -> {info['bam']}"
        )
    bam.close()
    params_by_level = {params.input_level: dataclasses.asdict(params)}
    paths = write_report(rows, args.out, sample, params_by_level)
    if args.emit_reads:
        reads_path = Path(args.out) / f"{sample}.kforce.reads.tsv"
        write_reads(read_rows, reads_path, sample)
        paths["reads"] = str(reads_path)
    if args.baseline_bam and args.emit_baseline_detail:
        info = write_baseline_detail(details, args.out, sample)
        baseline_bam.close()
        print(
            f"[kforce] patient reference: {info['n_diffs']} germline corrections, "
            f"{info['n_spliced']} observed-insert junctions -> {info['txt']}"
        )
    elif baseline_bam is not None:
        baseline_bam.close()

    n_pos = sum(1 for c in calls if c.call == "POSITIVE")
    n_unres = sum(1 for c in calls if c.call == "UNRESOLVED")
    unres_detail = f", {n_unres} UNRESOLVED" if n_unres else ""
    skip_detail = (
        (" | skipped " + ", ".join(f"{n} {kind}" for kind, n in skipped.most_common()))
        if skipped
        else ""
    )
    print(
        f"[kforce] {sample}: {len(calls)} variants force-called, {n_pos} POSITIVE"
        f"{unres_detail}{skip_detail} (ref={reference}, level={params.input_level}, "
        f"k={params.k}, D={params.max_mismatch})"
    )
    for c in calls:
        flag = "*" if c.call == "POSITIVE" else " "
        print(
            f"  {flag} {c.variant_id:42s} mut={c.n_mut:<4d} wt={c.n_wt:<4d} "
            f"vaf={c.vaf:.4f}  {c.call}"
        )
    print(f"[kforce] wrote: {paths['tidy']}")
    return 0


def cmd_inspect(args) -> int:
    """Run the 'inspect' subcommand: dump reconstructed sequences and k-mers."""
    params = _build_params(args)
    variants = load_manifest(args.manifest, patient=args.patient)
    if not variants:
        print("no variants loaded from manifest", file=sys.stderr)
        return 2
    name = args.name or Path(args.manifest).stem
    fa = pysam.FastaFile(args.reference)
    info = write_inspection(variants, fa, params, args.out, name)
    print(
        f"[kforce] inspected {info['n']} variants "
        f"({info['n_supported']} force-callable) -> {info['txt']}"
    )
    print(f"[kforce] machine-readable detail -> {info['tsv']}")
    return 0


def cmd_recall(args) -> int:
    """Run the 'recall' subcommand: re-derive calls with new thresholds."""
    from . import recall

    if Path(args.inp).is_dir():
        nf, nc = recall.recall_dir(args.inp, args.out, args.min_mut_molecules, args.vaf_floor)
        print(
            f"[recall] {nf} files, {nc} calls changed "
            f"(min_mut_molecules={args.min_mut_molecules}, vaf_floor={args.vaf_floor}) "
            f"-> {args.out}"
        )
    else:
        Path(args.out).mkdir(parents=True, exist_ok=True)
        out = str(Path(args.out) / Path(args.inp).name)
        nc = recall.recall_tsv(args.inp, out, args.min_mut_molecules, args.vaf_floor)
        print(f"[recall] {nc} calls changed -> {out}")
    return 0


def cmd_aggregate(args) -> int:
    """Run the 'aggregate' subcommand: concatenate per-sample tables."""
    inputs: list[str] = []
    for pat in args.inputs:
        inputs.extend(sorted(glob.glob(pat)) or [pat])
    out = Path(args.out)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=TIDY_COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        n = 0
        for path in inputs:
            with open(path) as ih:
                for row in csv.DictReader(ih, delimiter="\t"):
                    w.writerow(row)
                    n += 1
    print(f"[kforce] aggregated {len(inputs)} files, {n} rows -> {out}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser with all subcommands."""
    p = argparse.ArgumentParser(prog="kforce", description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)

    c = sub.add_parser("call", help="force-call manifest variants in one BAM/CRAM")
    c.add_argument("--bam", required=True, help="aligned, sorted BAM or CRAM")
    c.add_argument("--manifest", required=True, help="variant manifest (.xlsx)")
    c.add_argument("--reference", required=True, help="reference FASTA (indexed)")
    c.add_argument("--sample", help="sample name (default: BAM basename)")
    c.add_argument("--patient", help="filter manifest to this Patient")
    c.add_argument(
        "--level",
        default="wgs",
        choices=sorted(PRESETS),
        help="parameter preset / input_level label (default: wgs)",
    )
    c.add_argument("--out", default=".", help="output directory (default: .)")
    c.add_argument("--k", type=int)
    c.add_argument("--flank", type=int)
    c.add_argument("--window", type=int)
    c.add_argument("--max-mismatch", dest="max_mismatch", type=int)
    c.add_argument("--min-mut-kmers", dest="min_mut_kmers", type=int)
    c.add_argument("--min-wt-kmers", dest="min_wt_kmers", type=int)
    c.add_argument("--min-mut-molecules", dest="min_mut_molecules", type=int)
    c.add_argument("--vaf-floor", dest="vaf_floor", type=float)
    c.add_argument(
        "--min-indel-size",
        dest="min_indel_size",
        type=int,
        default=3,
        help="skip insertions/deletions whose size is <= this many bp (default 3: "
        "analyze only indels larger than 3 bp; 0 = load all). Symbolic <INS> "
        "(unknown size) and structural <DEL>/SV junctions are always kept.",
    )
    c.add_argument(
        "--max-reads-per-locus",
        dest="max_reads_per_locus",
        type=int,
        help="cap reads scanned per locus (0=unlimited); guards deep regions",
    )
    c.add_argument(
        "--min-mapq",
        dest="min_mapq",
        type=int,
        help="skip reads below this MAPQ when fetching fragments/baseline "
        "consensus (default 20; multi-mappers in repeats/segmental "
        "duplications produce spurious k-mer matches otherwise)",
    )
    c.add_argument(
        "--baseline-haplotype-margin",
        dest="baseline_haplotype_margin",
        type=int,
        help="+/- bp around a local indel's declared edit that --baseline-bam's "
        "observed_local_haplotype() compares baseline reads over (default "
        "20); empirically calibrated, see kforce/config.py",
    )
    c.add_argument(
        "--threads",
        type=int,
        default=1,
        help="htslib (de)compression threads for big BAM/CRAM reading",
    )
    c.add_argument(
        "--emit-reads",
        dest="emit_reads",
        action="store_true",
        help="also write <sample>.kforce.reads.tsv listing the MUT/WT "
        "supporting reads with their evidence",
    )
    c.add_argument(
        "--baseline-bam",
        dest="baseline_bam",
        help="patient baseline (diagnosis) BAM/CRAM; build patient-specific "
        "WT/MUT sequences from its pileup instead of the genome",
    )
    c.add_argument(
        "--baseline-reference",
        dest="baseline_reference",
        help="reference FASTA for the baseline BAM/CRAM when it was aligned to a "
        "different genome than --bam (e.g. targeted reads on one build, diagnosis "
        "CRAM on another). Defaults to --reference.",
    )
    c.add_argument(
        "--combined",
        action="store_true",
        help="two-tier (requires --baseline-bam): genome first, baseline only "
        "as a second-line check on genome-negatives; final call = genome "
        "if POSITIVE else baseline. Loss-free and compute-efficient.",
    )
    c.add_argument(
        "--emit-baseline-detail",
        dest="emit_baseline_detail",
        action="store_true",
        help="with --baseline-bam, write the patient "
        "reference review artifact (germline diffs, observed inserts)",
    )
    c.add_argument(
        "--emit-bam",
        dest="emit_bam",
        action="store_true",
        help="write <sample>.kforce.support.bam: MUT-supporting reads tagged "
        "KF:Z:<KMER|ANCHOR|BOTHEND|INS> plus ~10%% as many WT reads "
        "(KF:Z:WT) per variant, for IGV review (color/group by tag KF)",
    )
    c.add_argument(
        "--discordant-bothend",
        dest="bothend_parallel",
        action="store_true",
        help="for far-apart breakends, also count discordant read pairs (mates "
        "straddling the junction, no junction-spanning read) as MUT via the "
        "both-end signal. Off by default; raises chimera false-positive risk, "
        "so consider raising --min-mut-molecules with it.",
    )
    c.set_defaults(func=cmd_call)

    i = sub.add_parser(
        "inspect", help="dump per-variant detail (MUT/WT seqs, k-mers) for review; no BAM needed"
    )
    i.add_argument("--manifest", required=True, help="variant manifest (.xlsx)")
    i.add_argument("--reference", required=True, help="reference FASTA (indexed)")
    i.add_argument("--patient", help="filter manifest to this Patient")
    i.add_argument("--name", help="output basename (default: manifest stem)")
    i.add_argument(
        "--level",
        default="wgs",
        choices=sorted(PRESETS),
        help="parameter preset (controls k / flank / max_mismatch)",
    )
    i.add_argument("--out", default=".", help="output directory (default: .)")
    i.add_argument("--k", type=int)
    i.add_argument("--flank", type=int)
    i.add_argument("--max-mismatch", dest="max_mismatch", type=int)
    i.set_defaults(func=cmd_inspect)

    rc = sub.add_parser(
        "recall",
        help="re-derive calls from existing kforce.tsv with new "
        "min-mut-molecules / vaf-floor (no re-analysis)",
    )
    rc.add_argument("--in", dest="inp", required=True, help="kforce.tsv or directory of them")
    rc.add_argument("--out", required=True, help="output directory for recalled TSVs")
    rc.add_argument("--min-mut-molecules", dest="min_mut_molecules", type=int, required=True)
    rc.add_argument("--vaf-floor", dest="vaf_floor", type=float, default=0.0)
    rc.set_defaults(func=cmd_recall)

    a = sub.add_parser("aggregate", help="concatenate *.kforce.tsv into a cohort table")
    a.add_argument("inputs", nargs="+", help="tidy TSVs or globs")
    a.add_argument("--out", required=True, help="output cohort TSV")
    a.set_defaults(func=cmd_aggregate)
    return p


def main(argv=None) -> int:
    """Parse arguments and run the selected subcommand."""
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
