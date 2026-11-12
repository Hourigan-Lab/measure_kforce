"""End-to-end force-calling of manifest variants in one alignment file."""
from __future__ import annotations
import collections
import dataclasses
import pysam
from .alleles import SUPPORTED_KINDS, _anchors_far_apart, build_allele
from .baseline import ConsensusSource, build_patient_model
from .bothend import (
    JUNCTION_KINDS,
    classify_fragments_bothend,
    classify_fragments_parallel,
    is_bothend_candidate,
)
from .call import VariantCall, call_variant
from .classify import FragmentResult, classify_fragments
from .config import Params
from .manifest import Variant
from .reads import fetch_fragments, variant_loci


def _evidence_category(evidence: str) -> str:
    """Map an evidence path to the KF read-tag category."""
    return {"bothend": "BOTHEND", "novel_ins": "INS"}.get(evidence, "KMER")


def _load_collision_set(params: Params):
    """Load the optional genome-wide k-mer collision set."""
    if not params.genome_collision_path:
        return None
    import pickle

    with open(params.genome_collision_path, "rb") as f:
        return pickle.load(f)


def _fetches_distant_loci(v, params) -> bool:
    """Return True if a variant's reads come from loci too far apart to share one window."""
    loci = [(v.chrom, v.pos)]
    if v.breakend is not None:
        loci.append((v.breakend.mate_chrom, v.breakend.mate_pos))
    elif v.stop != v.pos:
        loci.append((v.chrom, v.stop))
    if len(loci) < 2:
        return False
    (c1, p1), (c2, p2) = loci[0], loci[1]
    return c1 != c2 or abs(p1 - p2) > params.window


def _bothend_safe_separation(v, params) -> bool:
    """Return True if the breakpoints are far enough apart for the both-end path."""
    if v.breakend is not None:
        c1, p1 = v.chrom, v.pos
        c2, p2 = v.breakend.mate_chrom, v.breakend.mate_pos
        if c1 != c2:
            return True
        sep = abs(p1 - p2)
    else:
        sep = abs(v.stop - v.pos)
    min_safe = 2 * (params.bothend_window + params.window)
    return sep > min_safe


def _has_long_homopolymer(km: str, max_run: int = 8) -> bool:
    """Return True if a k-mer contains a homopolymer run longer than a limit."""
    run = 1
    for i in range(1, len(km)):
        if km[i] == km[i - 1]:
            run += 1
            if run > max_run:
                return True
        else:
            run = 1
    return False


def _filter_model_kmers(model, genome_collision, params):
    """Remove genome-colliding and low-complexity k-mers from a model's MUT sets."""
    if not _fetches_distant_loci(model.variant, params):
        return model
    changed = False
    filtered_mut_unique = model.mut_unique
    filtered_fwd = model.mut_unique_fwd
    filtered_n = model.n_mut_unique
    if model.mut_unique:
        removed = set()
        if genome_collision:
            removed |= genome_collision
        removed |= {km for km in model.mut_unique if _has_long_homopolymer(km)}
        if removed:
            filtered_mut_unique = model.mut_unique - removed
            filtered_fwd = tuple(km for km in model.mut_unique_fwd if km not in removed)
            filtered_n = len(filtered_fwd)
            if len(filtered_mut_unique) < len(model.mut_unique):
                changed = True
    filtered_junction = model.mut_junction
    if model.mut_junction:
        removed_j = set()
        if genome_collision:
            removed_j |= genome_collision
        removed_j |= {km for km in model.mut_junction if _has_long_homopolymer(km)}
        if removed_j:
            filtered_junction = model.mut_junction - removed_j
            if len(filtered_junction) < len(model.mut_junction):
                changed = True
    if not changed:
        return model
    return dataclasses.replace(
        model,
        mut_unique=frozenset(filtered_mut_unique),
        mut_unique_fwd=filtered_fwd,
        n_mut_unique=filtered_n,
        mut_junction=frozenset(filtered_junction),
    )


def _try_delins_alt_trim(v, model, fragments, fa, params, genome_collision):
    """Retry a delins variant with single-base-trimmed ALTs when no MUT reads are found."""
    trims: list[str] = []
    if len(v.alt) > 1:
        trims.append(v.alt[1:])
        trims.append(v.alt[:-1])
    if len(v.alt) > 2:
        trims.append(v.alt[1:-1])

    best_model = model
    best_results = None
    best_n_mut = 0
    for trimmed_alt in trims:
        if not trimmed_alt:
            continue
        v_trimmed = dataclasses.replace(v, alt=trimmed_alt)
        try:
            model_trimmed = build_allele(v_trimmed, fa, params)
        except Exception:
            continue
        model_trimmed = _filter_model_kmers(model_trimmed, genome_collision, params)
        if model_trimmed.n_mut_unique < params.min_mut_kmers:
            continue
        results_trimmed = classify_fragments(fragments, model_trimmed, params)
        n_mut = sum(1 for r in results_trimmed.values() if r.label == "mut")
        if n_mut > best_n_mut:
            best_n_mut = n_mut
            best_model = model_trimmed
            best_results = results_trimmed

    if best_results is not None and best_n_mut >= params.min_mut_molecules:
        return best_model, best_results
    return model, None


def _detect_novel_ins_softclip(bam, v, fa, params):
    """Detect a long symbolic insertion from soft-clipped reads at the junction."""
    pos = v.pos
    chrom = v.chrom

    try:
        ref_downstream = fa.fetch(chrom, pos - 1, pos - 1 + 300).upper()
    except (ValueError, OSError):
        return set()

    supporting_qnames = set()
    window = 10
    start = max(0, pos - 1 - window)
    end = pos - 1 + window

    try:
        fetch_iter = bam.fetch(chrom, start, end)
    except (ValueError, OSError):
        return set()

    for r in fetch_iter:
        if r.is_secondary or r.is_supplementary:
            continue
        if r.mapping_quality < params.min_mapq:
            continue
        if r.cigartuples is None:
            continue

        seq = r.query_sequence
        if not seq:
            continue

        ref_pos = r.reference_start
        query_pos = 0

        for op, length in r.cigartuples:
            at_junction = abs(ref_pos - (pos - 1)) <= 5

            if op == 4 and at_junction and length >= 15:
                clip_seq = seq[query_pos : query_pos + length].upper()
                if len(clip_seq) >= 15:
                    check_seq = clip_seq[: min(20, len(clip_seq))]
                    if check_seq in ref_downstream:
                        supporting_qnames.add(r.query_name)

            if op in (0, 2, 3):
                ref_pos += length
            if op in (0, 1, 4):
                query_pos += length

    return supporting_qnames


def force_call(
    variants: list[Variant],
    fa: pysam.FastaFile,
    bam: pysam.AlignmentFile,
    params: Params,
    emit_reads: bool = False,
    model_builder=None,
    seq_source=None,
    emit_bam: bool = False,
) -> tuple[list[VariantCall], collections.Counter, list[dict], list, list]:
    """Reconstruct, fetch, classify and call every supported variant."""
    if model_builder is None:

        def model_builder(v):
            """Build the genome-based model for one variant."""
            return (build_allele(v, fa, params), None)

    if seq_source is None:
        seq_source = fa

    genome_collision = _load_collision_set(params)

    calls: list[VariantCall] = []
    skipped: collections.Counter = collections.Counter()
    read_rows: list[dict] = []
    details: list = []
    support: list = []
    for v in variants:
        if v.kind not in SUPPORTED_KINDS:
            skipped[v.kind] += 1
            continue
        model, detail = model_builder(v)
        if detail is not None:
            details.append(detail)
        model = _filter_model_kmers(model, genome_collision, params)
        fragments = fetch_fragments(bam, variant_loci(v), params)

        if model.novel_insertion:
            results = classify_fragments(fragments, model, params)
            evidence = "novel_ins"
            if not any(r.label == "mut" for r in results.values()):
                softclip_qnames = _detect_novel_ins_softclip(bam, v, fa, params)
                if len(softclip_qnames) >= params.min_mut_molecules:
                    for qn in softclip_qnames:
                        results[qn] = FragmentResult(
                            label="mut",
                            mut_votes=1,
                            wt_hits=0,
                            anchor_bridge=True,
                            anchor_gap=None,
                            mut_category="SOFTCLIP",
                        )
                    evidence = "softclip"
        elif is_bothend_candidate(v, model, params):
            results, _na, _nb = classify_fragments_bothend(
                fragments, v, model, seq_source, params, genome_collision
            )
            evidence = "bothend"
        elif (
            v.kind == "tandem_dup"
            and model.n_mut_unique == 0
            and model.head_anchor
            and model.tail_anchor
        ):
            results = classify_fragments(fragments, model, params)
            evidence = "anchor"
        elif (
            v.kind == "deletion"
            and model.n_mut_unique > 0
            and _anchors_far_apart(v, v.kind, params)
        ):
            model = dataclasses.replace(model, mut_junction=frozenset())
            results, _na, _nb = classify_fragments_parallel(
                fragments, v, model, seq_source, params, genome_collision
            )
            evidence = "kmer+bothend"
        elif (
            params.bothend_parallel
            and v.kind in JUNCTION_KINDS
            and model.n_mut_unique > 0
            and _anchors_far_apart(v, v.kind, params)
            and _bothend_safe_separation(v, params)
        ):
            results, _na, _nb = classify_fragments_parallel(
                fragments, v, model, seq_source, params, genome_collision
            )
            evidence = "kmer+bothend"
        else:
            results = classify_fragments(fragments, model, params)
            evidence = "kmer"
            if (
                v.kind == "delins"
                and model.n_mut_unique > 0
                and not any(r.label == "mut" for r in results.values())
            ):
                trimmed_model, trimmed_results = _try_delins_alt_trim(
                    v, model, fragments, fa, params, genome_collision
                )
                if trimmed_results is not None:
                    model = trimmed_model
                    results = trimmed_results
        calls.append(call_variant(model, results, params, evidence=evidence))

        if emit_reads:
            for qname, r in results.items():
                if r.label not in ("mut", "wt"):
                    continue
                row = {
                    "variant_id": v.variant_id,
                    "support": r.label.upper(),
                    "read_name": qname,
                    "evidence": evidence,
                    "wt_kmer_hits": r.wt_hits,
                    "read_seqs": ";".join(fragments[qname]),
                }
                if evidence == "bothend":
                    row.update(a_spec_hits=r.a_hits, b_spec_hits=r.b_hits)
                else:
                    row.update(
                        mut_kmer_votes=r.mut_votes,
                        anchor_bridge="yes" if r.anchor_bridge else "no",
                        anchor_gap="" if r.anchor_gap is None else r.anchor_gap,
                    )
                read_rows.append(row)

        if emit_bam:
            mut = {
                qn: (getattr(r, "mut_category", "") or _evidence_category(evidence))
                for qn, r in results.items()
                if r.label == "mut"
            }
            wt = [qn for qn, r in results.items() if r.label == "wt"]
            support.append({"vid": v.variant_id, "loci": variant_loci(v), "mut": mut, "wt": wt})
    return calls, skipped, read_rows, details, support


def force_call_combined(
    variants: list[Variant],
    fa: pysam.FastaFile,
    bam: pysam.AlignmentFile,
    params: Params,
    baseline_bam: pysam.AlignmentFile,
    emit_reads: bool = False,
    emit_bam: bool = False,
    baseline_fa: pysam.FastaFile = None,
):
    """Call with the genome first and the patient baseline only on genome-negatives."""
    g_calls, skipped, g_reads, _, g_support = force_call(
        variants, fa, bam, params, emit_reads=emit_reads, emit_bam=emit_bam
    )
    g_by = {c.variant_id: c for c in g_calls}

    neg = [
        v
        for v in variants
        if v.kind in SUPPORTED_KINDS
        and v.variant_id in g_by
        and g_by[v.variant_id].call != "POSITIVE"
    ]
    bfa = baseline_fa if baseline_fa is not None else fa
    ss = ConsensusSource(baseline_bam, bfa, params)

    def mb(v):
        """Build the patient-specific model for one variant."""
        return build_patient_model(v, baseline_bam, bfa, params)

    b_calls, _sk, b_reads, b_details, b_support = force_call(
        neg,
        fa,
        bam,
        params,
        emit_reads=emit_reads,
        model_builder=mb,
        seq_source=ss,
        emit_bam=emit_bam,
    )
    b_by = {c.variant_id: c for c in b_calls}

    calls, refs = [], []
    g_final, b_final = set(), set()
    for c in g_calls:
        if c.call != "POSITIVE" and c.variant_id in b_by:
            calls.append(b_by[c.variant_id])
            refs.append("baseline")
            b_final.add(c.variant_id)
        else:
            calls.append(c)
            refs.append("genome")
            g_final.add(c.variant_id)

    read_rows = []
    if emit_reads:
        read_rows = [r for r in g_reads if r["variant_id"] in g_final]
        read_rows += [r for r in b_reads if r["variant_id"] in b_final]
    support = []
    if emit_bam:
        support = [s for s in g_support if s["vid"] in g_final]
        support += [s for s in b_support if s["vid"] in b_final]
    return calls, refs, skipped, read_rows, b_details, support
