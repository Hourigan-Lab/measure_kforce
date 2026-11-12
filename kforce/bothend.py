"""Alignment-free both-end-flank detection for junctions inside repeats."""
from __future__ import annotations
from dataclasses import dataclass
import pysam
from .alleles import AlleleModel, _anchors_far_apart, _both_orient_kmers, _fetch
from .classify import FragmentResult, classify_fragments
from .config import Params
from .manifest import Variant
from .sequtil import is_acgt, iter_kmers


JUNCTION_KINDS = frozenset({"breakend", "tandem_dup", "deletion"})


def is_bothend_candidate(v: Variant, model: AlleleModel, params: Params) -> bool:
    """Return True if a variant should be resolved by the both-end fallback."""
    return (
        model.n_mut_unique == 0
        and v.kind in JUNCTION_KINDS
        and _anchors_far_apart(v, v.kind, params)
    )


def _breakpoint_loci(v: Variant) -> tuple[tuple[str, int], tuple[str, int]]:
    """Return the two breakpoint (chrom, pos) loci of a variant."""
    if v.breakend is not None:
        return (v.chrom, v.pos), (v.breakend.mate_chrom, v.breakend.mate_pos)
    return (v.chrom, v.pos), (v.chrom, v.stop)


def end_specific_kmers(
    v: Variant,
    fa: pysam.FastaFile,
    params: Params,
    genome_collision: frozenset[str] | None = None,
):
    """Return the k-mer sets unique to each end's flank."""
    (ca, pa), (cb, pb) = _breakpoint_loci(v)
    W, k = params.bothend_window, params.k
    ka = _both_orient_kmers(_fetch(fa, ca, pa - 1 - W, pa - 1 + W), k)
    kb = _both_orient_kmers(_fetch(fa, cb, pb - 1 - W, pb - 1 + W), k)
    a_spec = ka - kb
    b_spec = kb - ka
    if genome_collision is not None:
        a_spec = a_spec - genome_collision
        b_spec = b_spec - genome_collision
    return frozenset(a_spec), frozenset(b_spec)


@dataclass
class BothEndResult:
    """Classification of one fragment by the both-end method."""
    label: str
    a_hits: int
    b_hits: int
    wt_hits: int
    wt_hits_b: int = 0

    @property
    def mut_votes(self) -> int:
        """Return the weaker of the two end hit counts as the MUT vote."""
        return min(self.a_hits, self.b_hits)


def classify_fragment_bothend(
    seqs, a_spec, b_spec, model: AlleleModel, params: Params
) -> BothEndResult:
    """Label one fragment MUT/WT/ambiguous from its end-specific k-mer hits."""
    kms: set[str] = set()
    for s in seqs:
        for km in iter_kmers(s, params.k):
            if is_acgt(km):
                kms.add(km)
    a_hits = len(kms & a_spec)
    b_hits = len(kms & b_spec)
    wt_hits = len(kms & model.wt_span)
    wt_hits_b = len(kms & model.wt_span_b)
    t = params.bothend_min_kmers
    if a_hits >= t and b_hits >= t:
        label = "mut"
    elif max(wt_hits, wt_hits_b) >= params.min_wt_kmers:
        label = "wt"
    else:
        label = "ambiguous"
    return BothEndResult(
        label=label, a_hits=a_hits, b_hits=b_hits, wt_hits=wt_hits, wt_hits_b=wt_hits_b
    )


def classify_fragments_bothend(fragments, v, model, fa, params, genome_collision=None):
    """Classify all fragments with the both-end method."""
    a_spec, b_spec = end_specific_kmers(v, fa, params, genome_collision)
    results = {
        qn: classify_fragment_bothend(seqs, a_spec, b_spec, model, params)
        for qn, seqs in fragments.items()
    }
    return results, len(a_spec), len(b_spec)


def classify_fragments_parallel(fragments, v, model, fa, params, genome_collision=None):
    """Combine standard junction-k-mer evidence with both-end evidence per fragment."""
    kmer_res = classify_fragments(fragments, model, params)
    be_res, n_a, n_b = classify_fragments_bothend(fragments, v, model, fa, params, genome_collision)
    merged: dict[str, FragmentResult] = {}
    for qn, kr in kmer_res.items():
        br = be_res.get(qn)
        be_mut = br is not None and br.label == "mut"
        be_wt = br is not None and br.label == "wt"
        mut_category = ""
        if kr.label == "mut":
            label = "mut"
            mut_category = kr.mut_category
        elif be_mut:
            label = "mut"
            mut_category = "BOTHEND"
        elif kr.label == "wt" or be_wt:
            label = "wt"
        else:
            label = "ambiguous"
        merged[qn] = FragmentResult(
            label=label,
            mut_votes=max(kr.mut_votes, br.mut_votes if br else 0),
            wt_hits=max(kr.wt_hits, br.wt_hits if br else 0),
            anchor_bridge=kr.anchor_bridge,
            anchor_gap=kr.anchor_gap,
            mut_category=mut_category,
            wt_hits_b=max(kr.wt_hits_b, br.wt_hits_b if br else 0),
        )
    return merged, n_a, n_b
