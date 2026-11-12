"""Aggregate fragment classifications into per-variant counts and the MRD call."""
from __future__ import annotations
from dataclasses import asdict, dataclass
from typing import Optional
from .alleles import AlleleModel
from .config import Params, decide_call


@dataclass
class VariantCall:
    """Per-variant counts, VAF and call."""
    variant_id: str
    chrom: str
    pos: int
    mate_chrom: str
    mate_pos: str
    variant_type: str
    category: str
    patient: str
    n_mut: int
    n_wt: int
    n_ambiguous: int
    depth: int
    vaf: float
    call: str
    evidence: str
    n_mut_kmer_hits: int
    n_mut_unique: int
    manifest_t_vaf: float | None
    warnings: str
    vaf_bkpt_a: Optional[float] = None
    vaf_bkpt_b: Optional[float] = None

    def as_row(self, params: Params, sample: str, reference: str = "genome") -> dict:
        """Return the call as a flat dict for the tidy output table."""
        d = asdict(self)
        d.update(
            sample=sample,
            reference=reference,
            input_level=params.input_level,
            caller="kforce",
            k=params.k,
            max_mismatch=params.max_mismatch,
            min_mut_molecules=params.min_mut_molecules,
            vaf_floor=params.vaf_floor,
        )
        return d


def call_variant(
    model: AlleleModel,
    results: dict,
    params: Params,
    evidence: str = "kmer",
) -> VariantCall:
    """Count MUT/WT fragments for a variant and make its POSITIVE/NEGATIVE/UNRESOLVED call."""
    n_mut = sum(1 for r in results.values() if r.label == "mut")
    n_wt = sum(1 for r in results.values() if r.label == "wt")
    n_amb = sum(1 for r in results.values() if r.label == "ambiguous")
    depth = n_mut + n_wt
    mut_kmer_hits = sum(r.mut_votes for r in results.values() if r.label == "mut")

    n_wt_a = sum(
        1 for r in results.values() if r.label == "wt" and r.wt_hits >= getattr(r, "wt_hits_b", 0)
    )
    n_wt_b = n_wt - n_wt_a
    vaf_bkpt_a = vaf_bkpt_b = None
    has_mate = model.kind == "breakend" and bool(model.wt_span_b)
    if has_mate:
        vaf = (2 * n_mut) / (2 * n_mut + n_wt) if (2 * n_mut + n_wt) > 0 else 0.0
        vaf_bkpt_a = round(n_mut / (n_mut + n_wt_a), 6) if (n_mut + n_wt_a) > 0 else None
        vaf_bkpt_b = round(n_mut / (n_mut + n_wt_b), 6) if (n_mut + n_wt_b) > 0 else None
    else:
        vaf = n_mut / depth if depth > 0 else 0.0

    call = decide_call(
        n_mut, vaf, model.n_mut_unique, evidence, params.min_mut_molecules, params.vaf_floor
    )
    v = model.variant
    b = v.breakend
    return VariantCall(
        variant_id=v.variant_id,
        chrom=v.chrom,
        pos=v.pos,
        mate_chrom=b.mate_chrom if b else "",
        mate_pos=str(b.mate_pos) if b else "",
        variant_type=v.variant_type,
        category=v.category,
        patient=v.patient,
        n_mut=n_mut,
        n_wt=n_wt,
        n_ambiguous=n_amb,
        depth=depth,
        vaf=round(vaf, 6),
        call=call,
        evidence=evidence,
        n_mut_kmer_hits=mut_kmer_hits,
        n_mut_unique=model.n_mut_unique,
        manifest_t_vaf=v.t_vaf,
        warnings=";".join(model.warnings),
        vaf_bkpt_a=vaf_bkpt_a,
        vaf_bkpt_b=vaf_bkpt_b,
    )
