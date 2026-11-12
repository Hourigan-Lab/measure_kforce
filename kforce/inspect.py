"""Dump the sequences and k-mers that will be screened for each variant, for review."""
from __future__ import annotations
import csv
import textwrap
from pathlib import Path
import pysam
from .alleles import SUPPORTED_KINDS, UnsupportedVariant, build_allele
from .config import Params
from .manifest import Variant


def inspect_variant(v: Variant, fa: pysam.FastaFile, params: Params) -> dict:
    """Return a detail dict for one variant, supported or not."""
    d = {
        "variant_id": v.variant_id,
        "patient": v.patient,
        "kind": v.kind,
        "variant_type": v.variant_type,
        "category": v.category,
        "chrom": v.chrom,
        "pos": v.pos,
        "stop": v.stop,
        "ref": v.ref,
        "alt": v.alt,
        "t_vaf": v.t_vaf,
        "supported": v.kind in SUPPORTED_KINDS,
        "skip_reason": "",
        "mut_seq": "",
        "mut_len": 0,
        "junction_marked": "",
        "segments": [],
        "ref_seqs": [],
        "wt_span_seq": "",
        "n_mut_unique": 0,
        "mut_unique_kmers": [],
        "warnings": [],
    }
    if v.breakend is not None:
        b = v.breakend
        d["breakend"] = {
            "mate_chrom": b.mate_chrom,
            "mate_pos": b.mate_pos,
            "orient": b.orient.name,
            "insert": b.insert,
            "ref_base": b.ref_base,
        }
    if not d["supported"]:
        d["skip_reason"] = f"{v.kind} not force-called by k-force"
        return d

    try:
        m = build_allele(v, fa, params)
    except UnsupportedVariant as e:
        d["supported"] = False
        d["skip_reason"] = str(e)
        return d

    d.update(
        mut_seq=m.mut_seq,
        mut_len=len(m.mut_seq),
        junction_marked=m.junction_marked,
        segments=[(lab, seq) for lab, seq in m.segments],
        ref_seqs=[(lab, seq) for lab, seq in m.ref_seqs],
        wt_span_seq=m.wt_span_seq,
        n_mut_unique=m.n_mut_unique,
        mut_unique_kmers=list(m.mut_unique_fwd),
        warnings=list(m.warnings),
    )
    return d


def _wrap(seq: str, width: int = 80, indent: str = "      ") -> str:
    """Wrap a sequence to a fixed width with an indent."""
    if not seq:
        return indent + "(empty)"
    return "\n".join(indent + line for line in textwrap.wrap(seq, width))


def render_text(details: list[dict], params: Params) -> str:
    """Render the inspection details as a human-readable report."""
    out: list[str] = []
    n_sup = sum(1 for d in details if d["supported"])
    out.append("# kforce variant inspection")
    out.append(f"# params: k={params.k} flank={params.flank} max_mismatch={params.max_mismatch}")
    out.append(
        f"# {len(details)} variants ({n_sup} force-callable, {len(details) - n_sup} skipped)\n"
    )

    for d in details:
        out.append("=" * 88)
        out.append(f"variant_id: {d['variant_id']}")
        out.append(
            f"  patient={d['patient']}  type={d['variant_type']}  "
            f"category={d['category']}  kind={d['kind']}"
        )
        out.append(
            f"  locus: {d['chrom']}:{d['pos']}  Stop={d['stop']}  "
            f"Ref={d['ref']!r}  Alt={d['alt']!r}  manifest_T.VAF={d['t_vaf']}"
        )
        if "breakend" in d:
            b = d["breakend"]
            out.append(
                f"  breakend: mate={b['mate_chrom']}:{b['mate_pos']}  "
                f"orient={b['orient']}  insert={b['insert']!r}  "
                f"ref_base={b['ref_base']!r}"
            )
        if not d["supported"]:
            out.append(f"  supported: NO  ({d['skip_reason']})")
            out.append("")
            continue

        out.append("  supported: yes\n")
        out.append(f"  MUT haplotype  (len={d['mut_len']}, '|' marks junctions):")
        out.append(_wrap(d["junction_marked"], indent="    "))
        out.append("  MUT segments:")
        for i, (lab, seq) in enumerate(d["segments"], 1):
            out.append(f"    [{i}] {lab}  ({len(seq)} bp)")
            out.append(_wrap(seq, indent="        "))
        out.append("  WT / reference window(s) used for MUT-unique subtraction:")
        for lab, seq in d["ref_seqs"]:
            out.append(f"    {lab}  ({len(seq)} bp)")
            out.append(_wrap(seq, indent="        "))
        out.append(
            f"  WT-spanning sequence (defines the WT denominator, len={len(d['wt_span_seq'])}):"
        )
        out.append(_wrap(d["wt_span_seq"], indent="    "))
        out.append(f"  MUT-unique k-mers (k={params.k}): {d['n_mut_unique']}")
        for km in d["mut_unique_kmers"]:
            out.append(f"    {km}")
        out.append(f"  warnings: {'; '.join(d['warnings']) if d['warnings'] else '(none)'}")
        out.append("")
    return "\n".join(out)


_TSV_COLUMNS = [
    "variant_id",
    "patient",
    "kind",
    "variant_type",
    "category",
    "chrom",
    "pos",
    "stop",
    "ref",
    "alt",
    "t_vaf",
    "supported",
    "skip_reason",
    "mut_len",
    "n_mut_unique",
    "mut_seq",
    "junction_marked",
    "wt_span_seq",
    "ref_seqs",
    "mut_unique_kmers",
    "warnings",
]


def write_inspection(
    variants: list[Variant], fa: pysam.FastaFile, params: Params, outdir, name: str
) -> dict:
    """Write the inspection text and TSV files."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    details = [inspect_variant(v, fa, params) for v in variants]

    txt = outdir / f"{name}.kforce.inspect.txt"
    txt.write_text(render_text(details, params))

    tsv = outdir / f"{name}.kforce.inspect.tsv"
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=_TSV_COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for d in details:
            row = dict(d)
            row["ref_seqs"] = " ".join(f"{lab}={seq}" for lab, seq in d["ref_seqs"])
            row["mut_unique_kmers"] = ";".join(d["mut_unique_kmers"])
            row["warnings"] = "; ".join(d["warnings"])
            w.writerow(row)

    return {
        "txt": str(txt),
        "tsv": str(tsv),
        "n": len(details),
        "n_supported": sum(1 for d in details if d["supported"]),
    }
