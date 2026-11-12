"""Build patient-specific WT/MUT references from a baseline alignment instead of the genome."""
from __future__ import annotations
import collections
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import pysam
from .alleles import Reconstruction, _fetch, build_model, ref_and_guard_seqs, reconstruct
from .classify import _anchor_starts
from .config import Params
from .manifest import Variant
from .reads import fetch_fragments, variant_loci
from .sequtil import revcomp


class ConsensusSource:
    """Drop-in for a FASTA reader that returns the patient's per-position consensus sequence."""
    def __init__(self, bam: pysam.AlignmentFile, fa: pysam.FastaFile, params: Params):
        """Store the baseline alignment, genome reference and parameters."""
        self.bam = bam
        self.fa = fa
        self.params = params
        self._cache: dict[tuple[str, int, int], str] = {}

    def get_reference_length(self, chrom: str) -> int:
        """Return the length of a contig."""
        return self.fa.get_reference_length(chrom)

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Return the genome sequence for an interval with pileup consensus bases substituted."""
        key = (chrom, start, end)
        if key in self._cache:
            return self._cache[key]
        genome = self.fa.fetch(chrom, start, end).upper()
        out = list(genome)
        md, mf = self.params.baseline_min_depth, self.params.baseline_min_fraction
        if chrom in self.bam.references:
            try:
                for col in self.bam.pileup(
                    chrom,
                    start,
                    end,
                    truncate=True,
                    min_base_quality=0,
                    min_mapping_quality=self.params.min_mapq,
                ):
                    i = col.reference_pos - start
                    if i < 0 or i >= len(out):
                        continue
                    cnt = collections.Counter(
                        b.upper()
                        for b in col.get_query_sequences(add_indels=False)
                        if b and b.upper() in "ACGT"
                    )
                    tot = sum(cnt.values())
                    if tot < md:
                        continue
                    base, n = cnt.most_common(1)[0]
                    if n / tot >= mf:
                        out[i] = base
            except (OSError, ValueError):
                out = list(genome)
        res = "".join(out)
        self._cache[key] = res
        return res


@dataclass
class BaselineDetail:
    """Per-variant record of germline corrections and observed junction inserts."""
    variant_id: str
    germline_diffs: list[tuple[str, str, str]]
    manifest_insert: str
    observed_insert: Optional[str]
    spliced: bool
    n_mut_unique: int
    observed_note: str = ""


def _bridge_gap_seq(seq: str, head: str, tail: str, params: Params) -> Optional[str]:
    """Return the sequence between the head anchor's end and the tail anchor's start in a read."""
    k = len(head)
    heads = _anchor_starts(seq, head, params.max_mismatch)
    if not heads:
        return None
    tails = _anchor_starts(seq, tail, params.max_mismatch)
    if not tails:
        return None
    best = None
    best_abs = None
    for hi in heads:
        he = hi + k
        for ti in tails:
            gap = ti - he
            if -params.junction_hom_max <= gap <= params.junction_ins_max:
                a = abs(gap)
                if best_abs is None or a < best_abs:
                    best_abs = a
                    best = seq[he:ti] if gap > 0 else ""
    return best


def observed_insert(seqs, model, params: Params) -> Optional[str]:
    """Return the consensus untemplated insert at a junction across supporting baseline reads."""
    if not (model.head_anchor and model.tail_anchor):
        return None
    gaps = []
    for s in seqs:
        g = _bridge_gap_seq(s, model.head_anchor, model.tail_anchor, params)
        if g is None:
            g = _bridge_gap_seq(revcomp(s), model.head_anchor, model.tail_anchor, params)
        if g is not None:
            gaps.append(g)
    if len(gaps) < params.baseline_min_reads:
        return None
    common, n = collections.Counter(gaps).most_common(1)[0]
    return common if n >= params.baseline_min_reads else None


LOCAL_INDEL_KINDS = frozenset({"insertion", "deletion", "mnv", "delins"})


def _extract_local_haplotype(
    read: pysam.AlignedSegment, ref_start0: int, ref_end0: int
) -> Optional[str]:
    """Return a read's own sequence over a reference interval, or None if it does not span it."""
    if not read.cigartuples or not read.query_sequence:
        return None
    if read.reference_start > ref_start0:
        return None
    if read.reference_end is None or read.reference_end < ref_end0:
        return None
    rp, qp = read.reference_start, 0
    seq = read.query_sequence.upper()
    out: list[str] = []
    for op, ln in read.cigartuples:
        if op in (0, 7, 8):
            for _ in range(ln):
                if ref_start0 <= rp < ref_end0:
                    out.append(seq[qp])
                rp += 1
                qp += 1
        elif op == 1:
            if ref_start0 <= rp < ref_end0:
                out.append(seq[qp : qp + ln])
            qp += ln
        elif op in (2, 3):
            rp += ln
        elif op == 4:
            qp += ln
    return "".join(out)


def observed_local_haplotype(
    v: Variant,
    source: "ConsensusSource",
    bam: pysam.AlignmentFile,
    params: Params,
    margin: int = None,
):
    """Return the baseline-observed sequence for a local indel, or None if unchanged."""
    if margin is None:
        margin = params.baseline_haplotype_margin
    edit_start0, edit_end0 = v.pos - 1, v.stop
    ref_start0, ref_end0 = edit_start0 - margin, edit_end0 + margin
    ref_seq = _fetch(source, v.chrom, ref_start0, ref_end0)
    ins = "" if (v.alt in ("", "-") or v.alt.startswith("<")) else v.alt.upper()
    manifest_alt_seq = (
        _fetch(source, v.chrom, ref_start0, edit_start0)
        + ins
        + _fetch(source, v.chrom, edit_end0, ref_end0)
    )

    haps: collections.Counter = collections.Counter()
    try:
        for r in bam.fetch(v.chrom, max(0, ref_start0), ref_end0):
            if r.is_secondary or r.is_supplementary:
                continue
            if r.mapping_quality < params.min_mapq:
                continue
            h = _extract_local_haplotype(r, ref_start0, ref_end0)
            if h is not None:
                haps[h] += 1
    except (OSError, ValueError):
        return None

    if haps.get(manifest_alt_seq, 0) >= params.baseline_min_reads:
        return None

    candidates = [(h, n) for h, n in haps.items() if h not in (manifest_alt_seq, ref_seq)]
    if not candidates:
        return None
    candidates.sort(key=lambda hn: -hn[1])
    top_seq, top_n = candidates[0]
    runner_up_n = candidates[1][1] if len(candidates) > 1 else 0
    if top_n < params.baseline_min_reads or top_n < 2 * runner_up_n:
        return None

    note = (
        f"observed haplotype over {v.chrom}:{ref_start0 + 1}-{ref_end0} "
        f"({top_n} reads) differs from the manifest's predicted sequence "
        f"(0 exact matches; runner-up alternative had {runner_up_n} reads)"
    )
    return ref_start0, ref_end0, top_seq, note


def _germline_diffs(recon: Reconstruction, source: ConsensusSource, params: Params):
    """List positions where the patient consensus differs from the genome."""
    seen: dict[str, tuple[str, str, str]] = {}
    L = params.flank
    for label, pseq in recon.ref_seqs:
        chrom, posstr = label.rsplit(":", 1)
        c0 = int(posstr) - 1
        start = max(0, c0 - L)
        genome = source.fa.fetch(chrom, start, c0 + L).upper()
        for i, (g, p) in enumerate(zip(genome, pseq)):
            if g != p and p in "ACGT" and g in "ACGT":
                key = f"{chrom}:{start + i + 1}"
                seen[key] = (key, g, p)
    return sorted(seen.values(), key=lambda d: (d[0].split(":")[0], int(d[0].split(":")[1])))


def build_patient_model(
    v: Variant, baseline_bam: pysam.AlignmentFile, fa: pysam.FastaFile, params: Params
):
    """Build a patient-specific allele model for a variant from baseline reads."""
    source = ConsensusSource(baseline_bam, fa, params)
    recon = reconstruct(v, source, params)
    model = build_model(v, recon, params)

    diffs = _germline_diffs(recon, source, params)
    manifest_ins = v.breakend.insert if v.breakend is not None else ""
    observed = None
    spliced = False
    note = ""

    if model.head_anchor and model.tail_anchor:
        frags = fetch_fragments(baseline_bam, variant_loci(v), params)
        seqs = [s for reads in frags.values() for s in reads]
        observed = observed_insert(seqs, model, params)
        if observed is not None and observed != manifest_ins:
            flanks = (recon.segments[0], recon.segments[-1])
            new_segs = [flanks[0]]
            if observed:
                new_segs.append(("observed_insert", observed))
            new_segs.append(flanks[1])
            model = build_model(
                v,
                Reconstruction(
                    recon.kind,
                    new_segs,
                    recon.ref_seqs,
                    recon.wt_span_seq,
                    list(recon.warns),
                    guard_seqs=recon.guard_seqs,
                    wt_span_b_seq=recon.wt_span_b_seq,
                ),
                params,
            )
            spliced = True
    elif v.kind in LOCAL_INDEL_KINDS:
        obs = observed_local_haplotype(v, source, baseline_bam, params)
        if obs is not None:
            ref_start0, ref_end0, observed_seq, obs_note = obs
            L = params.flank
            left = _fetch(source, v.chrom, ref_start0 - L, ref_start0)
            right = _fetch(source, v.chrom, ref_end0, ref_end0 + L)
            new_segs = [
                (f"{v.chrom}:{ref_start0} left", left),
                ("observed_haplotype", observed_seq),
                (f"{v.chrom}:{ref_end0 + 1} right", right),
            ]
            new_ref_windows = [(v.chrom, ref_start0), (v.chrom, max(ref_start0, ref_end0 - 1))]
            new_ref_seqs, new_guard_seqs = ref_and_guard_seqs(source, new_ref_windows, params)
            model = build_model(
                v,
                Reconstruction(
                    recon.kind,
                    new_segs,
                    new_ref_seqs,
                    recon.wt_span_seq,
                    list(recon.warns),
                    guard_seqs=new_guard_seqs,
                ),
                params,
            )
            spliced = True
            observed = observed_seq
            note = obs_note

    if spliced and model.n_mut_unique > 0:
        model.novel_insertion = False

    detail = BaselineDetail(
        variant_id=v.variant_id,
        germline_diffs=diffs,
        manifest_insert=manifest_ins,
        observed_insert=observed,
        spliced=spliced,
        n_mut_unique=model.n_mut_unique,
        observed_note=note,
    )
    return model, detail


def write_baseline_detail(details: list[BaselineDetail], outdir, sample: str) -> dict:
    """Write the patient-reference review artifacts (TSV and text)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tsv = outdir / f"{sample}.kforce.baseline_detail.tsv"
    txt = outdir / f"{sample}.kforce.baseline_detail.txt"
    cols = [
        "variant_id",
        "n_germline_diffs",
        "germline_diffs",
        "manifest_insert",
        "observed_insert",
        "spliced",
        "observed_note",
        "n_mut_unique",
    ]
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for d in details:
            w.writerow(
                {
                    "variant_id": d.variant_id,
                    "n_germline_diffs": len(d.germline_diffs),
                    "germline_diffs": ";".join(f"{p}:{g}>{q}" for p, g, q in d.germline_diffs),
                    "manifest_insert": d.manifest_insert,
                    "observed_insert": "" if d.observed_insert is None else d.observed_insert,
                    "spliced": "yes" if d.spliced else "no",
                    "observed_note": d.observed_note,
                    "n_mut_unique": d.n_mut_unique,
                }
            )
    n_diffs = sum(len(d.germline_diffs) for d in details)
    n_spliced = sum(1 for d in details if d.spliced)
    lines = [
        f"# kforce patient baseline reference — {sample}",
        f"# {len(details)} variants | {n_diffs} germline corrections | "
        f"{n_spliced} junctions with an observed insert\n",
    ]
    for d in details:
        lines.append(f"{d.variant_id}")
        gl = "; ".join(f"{p} {g}>{q}" for p, g, q in d.germline_diffs) or "(none)"
        lines.append(f"  germline diffs vs genome ({len(d.germline_diffs)}): {gl}")
        if d.manifest_insert or d.observed_insert is not None:
            obs = "(not observed)" if d.observed_insert is None else repr(d.observed_insert)
            tag = "  [SPLICED into MUT]" if d.spliced else ""
            lines.append(f"  junction insert: manifest={d.manifest_insert!r} observed={obs}{tag}")
        if d.observed_note:
            lines.append(f"  local-indel correction: {d.observed_note}  [SPLICED into MUT]")
        lines.append(f"  MUT-unique k-mers: {d.n_mut_unique}\n")
    txt.write_text("\n".join(lines))
    return {"tsv": str(tsv), "txt": str(txt), "n_diffs": n_diffs, "n_spliced": n_spliced}
