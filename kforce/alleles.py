"""Reconstruct local WT / MUT sequences for a variant and build its k-mer dictionaries."""
from __future__ import annotations
from dataclasses import dataclass
import pysam
from .config import Params
from .manifest import BndOrient, Variant
from .sequtil import is_acgt, iter_kmers, min_hamming, revcomp

SUPPORTED_KINDS = frozenset({"breakend", "tandem_dup", "insertion", "deletion", "mnv", "delins"})


class UnsupportedVariant(Exception):
    """Raised for variants that cannot be force-called (SNVs, unknown symbolic ALTs)."""


@dataclass
class AlleleModel:
    """Holds a variant's MUT/WT sequences, k-mer dictionaries and anchors."""
    variant: Variant
    mut_seq: str
    wt_span_seq: str
    mut_unique: frozenset[str]
    wt_span: frozenset[str]
    ref_all: frozenset[str]
    n_mut_unique: int
    mut_junction: frozenset[str] = frozenset()
    wt_span_b: frozenset[str] = frozenset()
    kind: str = ""
    segments: tuple[tuple[str, str], ...] = ()
    ref_seqs: tuple[tuple[str, str], ...] = ()
    mut_unique_fwd: tuple[str, ...] = ()
    head_anchor: str = ""
    tail_anchor: str = ""
    anchor_source: str = ""
    novel_insertion: bool = False
    warnings: tuple[str, ...] = ()

    @property
    def junction_marked(self) -> str:
        """Return the MUT sequence with '|' marking segment boundaries."""
        return "|".join(s for _, s in self.segments)


def _fetch(fa: pysam.FastaFile, chrom: str, start: int, end: int) -> str:
    """Fetch an upper-cased, contig-clamped 0-based half-open reference interval."""
    length = fa.get_reference_length(chrom)
    start = max(0, start)
    end = min(length, end)
    if start >= end:
        return ""
    return fa.fetch(chrom, start, end).upper()


def _both_orient_kmers(seq: str, k: int) -> set[str]:
    """Return all ACGT-only k-mers of a sequence plus their reverse complements."""
    out: set[str] = set()
    for km in iter_kmers(seq, k):
        if is_acgt(km):
            out.add(km)
            out.add(revcomp(km))
    return out


Segments = list


def _drop_empty(segs: Segments) -> Segments:
    """Remove empty segments from a labelled segment list."""
    return [(lab, s) for lab, s in segs if s]


def _reconstruct_bnd(v: Variant, fa: pysam.FastaFile, L: int) -> Segments:
    """Build the labelled MUT junction segments for a breakend variant."""
    b = v.breakend
    p0 = v.pos - 1
    m0 = b.mate_pos - 1
    ins = b.insert
    loc = f"{v.chrom}:{v.pos}"
    mate = f"{b.mate_chrom}:{b.mate_pos}"

    if b.orient is BndOrient.LEFT_FWD:
        return _drop_empty(
            [
                (f"{loc} left", _fetch(fa, v.chrom, p0 - L + 1, p0 + 1)),
                ("insert", ins),
                (f"{mate} mate(+)", _fetch(fa, b.mate_chrom, m0, m0 + L)),
            ]
        )
    if b.orient is BndOrient.LEFT_REV:
        return _drop_empty(
            [
                (f"{loc} left", _fetch(fa, v.chrom, p0 - L + 1, p0 + 1)),
                ("insert", ins),
                (
                    f"{mate} mate(rc)",
                    revcomp(_fetch(fa, b.mate_chrom, m0 - L + 1, m0 + 1)),
                ),
            ]
        )
    if b.orient is BndOrient.RIGHT_FWD:
        return _drop_empty(
            [
                (f"{mate} mate(+)", _fetch(fa, b.mate_chrom, m0 - L + 1, m0 + 1)),
                ("insert", ins),
                (f"{loc} right", _fetch(fa, v.chrom, p0, p0 + L)),
            ]
        )
    if b.orient is BndOrient.RIGHT_REV:
        return _drop_empty(
            [
                (f"{mate} mate(rc)", revcomp(_fetch(fa, b.mate_chrom, m0, m0 + L))),
                ("insert", ins),
                (f"{loc} right", _fetch(fa, v.chrom, p0, p0 + L)),
            ]
        )
    raise AssertionError(f"unhandled orient {b.orient}")


def _reconstruct_smallvar(v: Variant, fa: pysam.FastaFile, L: int) -> Segments:
    """Build the labelled MUT segments for an indel, delins or MNV."""
    ins = "" if (v.alt in ("", "-") or v.alt.startswith("<")) else v.alt.upper()
    ins_label = "insert(ALT)" if ins else "insert"
    return _drop_empty(
        [
            (
                f"{v.chrom}:{v.pos} left",
                _fetch(fa, v.chrom, (v.pos - 1) - L, v.pos - 1),
            ),
            (ins_label, ins),
            (f"{v.chrom}:{v.stop} right", _fetch(fa, v.chrom, v.stop, v.stop + L)),
        ]
    )


def _reconstruct_tandem_dup(v: Variant, fa: pysam.FastaFile, L: int) -> Segments:
    """Build the labelled head-to-tail MUT junction for a tandem duplication."""
    return _drop_empty(
        [
            (f"copy-end @{v.chrom}:{v.stop}", _fetch(fa, v.chrom, v.stop - L, v.stop)),
            (
                f"copy-start @{v.chrom}:{v.pos}",
                _fetch(fa, v.chrom, v.pos - 1, v.pos - 1 + L),
            ),
        ]
    )


def _anchors_far_apart(v: Variant, kind: str, params: Params) -> bool:
    """Return True if the two boundary anchors are far enough apart to be used safely."""
    sep = params.anchor_min_separation
    if kind == "breakend":
        b = v.breakend
        return b.mate_chrom != v.chrom or abs(v.pos - b.mate_pos) > sep
    if kind == "tandem_dup":
        return abs(v.stop - v.pos) > sep
    if kind == "deletion":
        return abs(v.stop - v.pos) > sep
    return False


def _anchor_pair_safe(v: Variant, kind: str, params: Params) -> bool:
    """Return True if the order-based anchor-pair bridge cannot be mimicked by a wild-type read."""
    if kind == "tandem_dup":
        return abs(v.stop - v.pos) >= 2 * params.k
    return _anchors_far_apart(v, kind, params)


@dataclass
class Reconstruction:
    """Sequences read from a source (genome or consensus) for one variant."""
    kind: str
    segments: list[tuple[str, str]]
    ref_seqs: list[tuple[str, str]]
    wt_span_seq: str
    warns: list[str]
    guard_seqs: tuple[tuple[str, str], ...] = ()
    wt_span_b_seq: str = ""


def ref_and_guard_seqs(
    source, ref_windows: list[tuple[str, int]], params: Params
) -> tuple[list[tuple[str, str]], tuple[tuple[str, str], ...]]:
    """Fetch the local and broad reference windows around each locus of a variant."""
    L = params.flank
    ref_seqs = [(f"{c}:{c0 + 1}", _fetch(source, c, c0 - L, c0 + L)) for c, c0 in ref_windows]
    W = max(params.window, L)
    guard_seqs = tuple(
        (f"{c}:{c0 + 1}", _fetch(source, c, c0 - W, c0 + W)) for c, c0 in ref_windows
    )
    return ref_seqs, guard_seqs


def reconstruct(v: Variant, source, params: Params) -> Reconstruction:
    """Reconstruct the MUT segments and reference windows for a variant from a sequence source."""
    kind = v.kind
    if kind not in SUPPORTED_KINDS:
        raise UnsupportedVariant(f"{v.variant_id}: kind={kind}")
    L = params.flank
    warns: list[str] = []

    if kind == "breakend":
        segments = _reconstruct_bnd(v, source, L)
        ref_windows = [
            (v.chrom, v.pos - 1),
            (v.breakend.mate_chrom, v.breakend.mate_pos - 1),
        ]
        wt_chrom, wt0 = v.chrom, v.pos - 1
        ref_at_pos = _fetch(source, v.chrom, v.pos - 1, v.pos)
        if v.ref and ref_at_pos and ref_at_pos != v.ref[0]:
            warns.append(
                f"reference base {ref_at_pos!r} at {v.chrom}:{v.pos} != manifest REF {v.ref!r}"
            )
    elif kind == "tandem_dup":
        segments = _reconstruct_tandem_dup(v, source, L)
        ref_windows = [(v.chrom, v.pos - 1), (v.chrom, v.stop - 1)]
        wt_chrom, wt0 = v.chrom, v.stop - 1
    else:
        segments = _reconstruct_smallvar(v, source, L)
        ref_windows = [(v.chrom, v.pos - 1)]
        if v.stop != v.pos:
            ref_windows.append((v.chrom, v.stop - 1))
        wt_chrom, wt0 = v.chrom, v.pos - 1

    ref_seqs, guard_seqs = ref_and_guard_seqs(source, ref_windows, params)
    wt_span_seq = _fetch(source, wt_chrom, wt0 - (params.k - 1), wt0 + params.k)
    wt_span_b_seq = ""
    if kind == "breakend":
        m0 = v.breakend.mate_pos - 1
        wt_span_b_seq = _fetch(source, v.breakend.mate_chrom, m0 - (params.k - 1), m0 + params.k)
    return Reconstruction(
        kind=kind,
        segments=segments,
        ref_seqs=ref_seqs,
        wt_span_seq=wt_span_seq,
        warns=warns,
        guard_seqs=guard_seqs,
        wt_span_b_seq=wt_span_b_seq,
    )


def build_model(v: Variant, recon: Reconstruction, params: Params) -> AlleleModel:
    """Build the k-mer dictionaries and anchors from a Reconstruction without any I/O."""
    k = params.k
    segments = recon.segments
    warns = list(recon.warns)
    mut_seq = "".join(s for _, s in segments)

    ref_kmers: set[str] = set()
    for _, seq in recon.ref_seqs:
        ref_kmers |= _both_orient_kmers(seq, k)
    ref_guard = set(ref_kmers)
    for _, seq in recon.guard_seqs:
        ref_guard |= _both_orient_kmers(seq, k)
    wt_span = _both_orient_kmers(recon.wt_span_seq, k)
    wt_span_b = _both_orient_kmers(recon.wt_span_b_seq, k) if recon.wt_span_b_seq else set()

    margin = 2 * params.max_mismatch
    mut_unique_fwd: set[str] = set()
    mut_junction_fwd: set[str] = set()
    for km in iter_kmers(mut_seq, k):
        if not is_acgt(km):
            continue
        if km in ref_kmers:
            continue
        if min_hamming(km, ref_kmers) <= margin:
            if km not in ref_guard:
                mut_junction_fwd.add(km)
            continue
        mut_unique_fwd.add(km)

    mut_unique = set(mut_unique_fwd)
    for km in mut_unique_fwd:
        mut_unique.add(revcomp(km))

    mut_junction = set(mut_junction_fwd)
    for km in mut_junction_fwd:
        mut_junction.add(revcomp(km))

    if not mut_unique_fwd:
        warns.append("no MUT-unique k-mers (low complexity / repeat region?)")

    head_anchor = tail_anchor = ""
    anchor_source = ""
    if (
        recon.kind == "tandem_dup"
        and len(segments) >= 2
        and _anchor_pair_safe(v, recon.kind, params)
    ):
        h = segments[0][1][-k:]
        t = segments[-1][1][:k]
        if len(h) == k and is_acgt(h):
            head_anchor = h
        if len(t) == k and is_acgt(t):
            tail_anchor = t
        anchor_source = "reference" if not mut_unique_fwd else "mut_unique"
    elif mut_unique_fwd and len(segments) >= 2 and _anchor_pair_safe(v, recon.kind, params):
        h = segments[0][1][-k:]
        t = segments[-1][1][:k]
        if len(h) == k and is_acgt(h):
            head_anchor = h
        if len(t) == k and is_acgt(t):
            tail_anchor = t
        anchor_source = "mut_unique"

    novel_insertion = recon.kind == "insertion" and v.alt.startswith("<") and len(segments) >= 2
    if novel_insertion:
        h = segments[0][1][-k:]
        t = segments[-1][1][:k]
        head_anchor = h if (len(h) == k and is_acgt(h)) else ""
        tail_anchor = t if (len(t) == k and is_acgt(t)) else ""
        novel_insertion = bool(head_anchor and tail_anchor)

    return AlleleModel(
        variant=v,
        mut_seq=mut_seq,
        wt_span_seq=recon.wt_span_seq,
        mut_unique=frozenset(mut_unique),
        wt_span=frozenset(wt_span),
        wt_span_b=frozenset(wt_span_b),
        ref_all=frozenset(ref_guard),
        n_mut_unique=len(mut_unique_fwd),
        mut_junction=frozenset(mut_junction)
        if len(mut_unique_fwd) < params.min_mut_kmers
        else frozenset(),
        kind=recon.kind,
        segments=tuple(segments),
        ref_seqs=tuple(recon.ref_seqs),
        mut_unique_fwd=tuple(sorted(mut_unique_fwd)),
        novel_insertion=novel_insertion,
        head_anchor=head_anchor,
        tail_anchor=tail_anchor,
        anchor_source=anchor_source,
        warnings=tuple(warns),
    )


def build_allele(v: Variant, fa: pysam.FastaFile, params: Params) -> AlleleModel:
    """Reconstruct sequences and build the k-mer model for one variant from a reference FASTA."""
    return build_model(v, reconstruct(v, fa, params), params)
