"""Classify each fragment as MUT-supporting, WT-supporting or ambiguous."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Optional
from .alleles import AlleleModel
from .config import Params
from .sequtil import hamming, is_acgt, iter_kmers, min_hamming, revcomp


@dataclass
class FragmentResult:
    """Classification and supporting evidence for one fragment."""
    label: str
    mut_votes: int
    wt_hits: int
    anchor_bridge: bool = False
    anchor_gap: Optional[int] = None
    mut_category: str = ""
    wt_hits_b: int = 0


def _anchor_starts(seq: str, anchor: str, max_mm: int) -> list[int]:
    """Return start indices where a window matches the anchor within a mismatch limit."""
    k = len(anchor)
    out = []
    for i in range(len(seq) - k + 1):
        if hamming(seq[i : i + k], anchor) <= max_mm:
            out.append(i)
    return out


def _bridge_gap(seq: str, head: str, tail: str, params: Params) -> Optional[int]:
    """Return the smallest gap at which the head anchor is followed by the tail anchor."""
    k = len(head)
    heads = _anchor_starts(seq, head, params.anchor_max_mismatch)
    if not heads:
        return None
    tails = _anchor_starts(seq, tail, params.anchor_max_mismatch)
    if not tails:
        return None
    best = None
    for hi in heads:
        head_end = hi + k
        for ti in tails:
            gap = ti - head_end
            if -params.junction_hom_max <= gap <= params.junction_ins_max:
                if best is None or abs(gap) < abs(best):
                    best = gap
    return best


def classify_fragment(seqs: list[str], model: AlleleModel, params: Params) -> FragmentResult:
    """Classify one fragment by MUT-unique k-mer votes, WT k-mers and the anchor-pair bridge."""
    k = params.k
    D = params.max_mismatch

    read_kmers: set[str] = set()
    for s in seqs:
        for km in iter_kmers(s, k):
            if is_acgt(km):
                read_kmers.add(km)

    mut_votes = 0
    wt_hits = 0
    wt_hits_b = 0
    for km in read_kmers:
        if km in model.mut_unique:
            mut_votes += 1
            continue
        if km in model.mut_junction:
            mut_votes += 1
            continue
        if km in model.wt_span:
            wt_hits += 1
            continue
        if km in model.wt_span_b:
            wt_hits_b += 1
            continue
        if km in model.ref_all:
            continue
        if D > 0 and min_hamming(km, model.mut_unique) <= D:
            mut_votes += 1

    anchor_gap = None
    if model.head_anchor and model.tail_anchor:
        for s in seqs:
            g = _bridge_gap(s, model.head_anchor, model.tail_anchor, params)
            if g is None:
                g = _bridge_gap(revcomp(s), model.head_anchor, model.tail_anchor, params)
            if g is not None:
                anchor_gap = g
                break
    anchor_bridge = anchor_gap is not None

    mut_category = ""
    if mut_votes >= params.min_mut_kmers:
        label = "mut"
        mut_category = "KMER"
    elif anchor_bridge:
        label = "mut"
        mut_category = "ANCHOR"
    elif max(wt_hits, wt_hits_b) >= params.min_wt_kmers:
        label = "wt"
    else:
        label = "ambiguous"
    return FragmentResult(
        label=label,
        mut_votes=mut_votes,
        wt_hits=wt_hits,
        anchor_bridge=anchor_bridge,
        anchor_gap=anchor_gap,
        mut_category=mut_category,
        wt_hits_b=wt_hits_b,
    )


def classify_fragment_novel_ins(seqs, model: AlleleModel, params: Params) -> FragmentResult:
    """Classify one fragment for a symbolic insertion by the gap between its two flank anchors."""
    head, tail, k, D = (
        model.head_anchor,
        model.tail_anchor,
        len(model.head_anchor),
        params.anchor_max_mismatch,
    )
    mut_gap = None
    wt = False
    for s in seqs:
        for seq in (s, revcomp(s)):
            heads = _anchor_starts(seq, head, D)
            if not heads:
                continue
            tails = _anchor_starts(seq, tail, D)
            if not tails:
                continue
            for hi in heads:
                he = hi + k
                for ti in tails:
                    gap = ti - he
                    if gap >= params.min_ins_gap:
                        mut_gap = gap if mut_gap is None else min(mut_gap, gap)
                    elif -params.junction_hom_max <= gap < params.min_ins_gap:
                        wt = True
    if mut_gap is not None:
        return FragmentResult(
            "mut",
            mut_votes=1,
            wt_hits=0,
            anchor_bridge=True,
            anchor_gap=mut_gap,
            mut_category="INS",
        )
    if wt:
        return FragmentResult("wt", mut_votes=0, wt_hits=1)
    return FragmentResult("ambiguous", mut_votes=0, wt_hits=0)


def classify_fragments(
    fragments: dict[str, list[str]], model: AlleleModel, params: Params
) -> dict[str, FragmentResult]:
    """Classify every fragment of a variant, caching identical read sets."""
    fn = classify_fragment_novel_ins if model.novel_insertion else classify_fragment
    cache: dict[tuple, FragmentResult] = {}
    out: dict[str, FragmentResult] = {}
    for qn, seqs in fragments.items():
        key = tuple(sorted(seqs))
        r = cache.get(key)
        if r is None:
            r = fn(seqs, model, params)
            cache[key] = r
        out[qn] = r
    return out