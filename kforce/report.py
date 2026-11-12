"""Write the per-sample TSV, summary, HTML and JSON outputs."""
from __future__ import annotations
import csv
import html
import json
from pathlib import Path

TIDY_COLUMNS = [
    "sample",
    "patient",
    "variant_id",
    "caller",
    "reference",
    "input_level",
    "chrom",
    "pos",
    "mate_chrom",
    "mate_pos",
    "variant_type",
    "category",
    "n_mut",
    "n_wt",
    "n_ambiguous",
    "depth",
    "vaf",
    "vaf_bkpt_a",
    "vaf_bkpt_b",
    "call",
    "evidence",
    "manifest_t_vaf",
    "n_mut_kmer_hits",
    "n_mut_unique",
    "k",
    "max_mismatch",
    "min_mut_molecules",
    "vaf_floor",
    "warnings",
]


def write_tidy(rows: list[dict], path: Path) -> None:
    """Write the tidy one-row-per-variant TSV."""
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=TIDY_COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _summary_table(rows: list[dict]) -> tuple[list[str], list[dict]]:
    """Pivot tidy rows to one record per variant with columns per input level."""
    levels = sorted({r["input_level"] for r in rows})
    by_var: dict[str, dict] = {}
    for r in rows:
        rec = by_var.setdefault(
            r["variant_id"],
            {
                "variant_id": r["variant_id"],
                "variant_type": r["variant_type"],
                "manifest_t_vaf": r.get("manifest_t_vaf"),
            },
        )
        lv = r["input_level"]
        rec[f"{lv}_mut"] = r["n_mut"]
        rec[f"{lv}_wt"] = r["n_wt"]
        rec[f"{lv}_vaf"] = r["vaf"]
        rec[f"{lv}_call"] = r["call"]
    cols = ["variant_id", "variant_type", "manifest_t_vaf"]
    for lv in levels:
        cols += [f"{lv}_mut", f"{lv}_wt", f"{lv}_vaf", f"{lv}_call"]
    return cols, list(by_var.values())


def write_summary_tsv(rows: list[dict], path: Path) -> None:
    """Write the wide per-variant summary TSV."""
    cols, recs = _summary_table(rows)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for rec in recs:
            w.writerow(rec)


def write_summary_html(rows: list[dict], path: Path, sample: str) -> None:
    """Write the per-variant summary as an HTML table."""
    cols, recs = _summary_table(rows)
    th = "".join(f"<th>{html.escape(c)}</th>" for c in cols)
    body = []
    for rec in recs:
        tds = []
        for c in cols:
            val = rec.get(c, "")
            cell = "" if val is None else html.escape(str(val))
            if c.endswith("_call"):
                color = "#c0392b" if val == "POSITIVE" else "#7f8c8d"
                tds.append(f'<td style="font-weight:600;color:{color}">{cell}</td>')
            else:
                tds.append(f"<td>{cell}</td>")
        body.append("<tr>" + "".join(tds) + "</tr>")
    doc = f"""<!doctype html><meta charset="utf-8">
<title>kforce — {html.escape(sample)}</title>
<style>
 body{{font-family:system-ui,Arial,sans-serif;margin:1.5rem}}
 table{{border-collapse:collapse;font-size:13px}}
 th,td{{border:1px solid #ddd;padding:4px 8px;text-align:right}}
 th{{background:#f4f6f8}} td:first-child,th:first-child{{text-align:left}}
 tr:nth-child(even){{background:#fafbfc}}
</style>
<h2>kforce force-calling — {html.escape(sample)}</h2>
<table><thead><tr>{th}</tr></thead><tbody>{"".join(body)}</tbody></table>
"""
    path.write_text(doc)


def write_json(rows: list[dict], path: Path, sample: str, params_by_level: dict) -> None:
    """Write nested sample/variant/level metrics and parameters as JSON."""
    variants: dict[str, dict] = {}
    for r in rows:
        v = variants.setdefault(
            r["variant_id"],
            {
                "chrom": r["chrom"],
                "pos": r["pos"],
                "mate_chrom": r["mate_chrom"],
                "mate_pos": r["mate_pos"],
                "variant_type": r["variant_type"],
                "manifest_t_vaf": r.get("manifest_t_vaf"),
                "levels": {},
            },
        )
        v["levels"][r["input_level"]] = {
            "n_mut": r["n_mut"],
            "n_wt": r["n_wt"],
            "n_ambiguous": r["n_ambiguous"],
            "depth": r["depth"],
            "vaf": r["vaf"],
            "call": r["call"],
            "n_mut_kmer_hits": r["n_mut_kmer_hits"],
            "n_mut_unique": r["n_mut_unique"],
        }
    doc = {"sample": sample, "params": params_by_level, "variants": variants}
    path.write_text(json.dumps(doc, indent=2))


READS_COLUMNS = [
    "sample",
    "variant_id",
    "support",
    "read_name",
    "evidence",
    "mut_kmer_votes",
    "anchor_bridge",
    "anchor_gap",
    "a_spec_hits",
    "b_spec_hits",
    "wt_kmer_hits",
    "read_seqs",
]


def write_reads(read_rows: list[dict], path: Path, sample: str) -> None:
    """Write per-read MUT/WT evidence."""
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=READS_COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in read_rows:
            w.writerow({**r, "sample": sample})


def write_report(rows: list[dict], outdir, sample: str, params_by_level: dict) -> dict:
    """Write all per-sample outputs and return their paths."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = {
        "tidy": outdir / f"{sample}.kforce.tsv",
        "summary": outdir / f"{sample}.kforce.summary.tsv",
        "html": outdir / f"{sample}.kforce.summary.html",
        "json": outdir / f"{sample}.kforce.json",
    }
    write_tidy(rows, paths["tidy"])
    write_summary_tsv(rows, paths["summary"])
    write_summary_html(rows, paths["html"], sample)
    write_json(rows, paths["json"], sample, params_by_level)
    return {k: str(v) for k, v in paths.items()}
