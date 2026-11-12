"""Run parameters (Params), presets and the shared call-decision rule."""
from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True)
class Params:
    """Immutable set of run parameters shared by every module."""
    k: int = 25
    flank: int = 150
    min_anchor: int = 5

    max_mismatch: int = 1

    window: int = 600
    max_reads_per_locus: int = 0
    min_mapq: int = 20

    min_mut_kmers: int = 3
    min_wt_kmers: int = 3

    junction_ins_max: int = 30
    junction_hom_max: int = 10
    anchor_max_mismatch: int | None = None
    min_ins_gap: int = 1
    anchor_min_separation: int = 300

    bothend_window: int = 400
    bothend_min_kmers: int = 30
    bothend_parallel: bool = False

    genome_collision_path: str | None = None

    baseline_min_depth: int = 5
    baseline_min_fraction: float = 0.6
    baseline_min_reads: int = 2
    baseline_haplotype_margin: int = 20

    min_mut_molecules: int = 2
    vaf_floor: float = 0.0

    input_level: str = "wgs"

    def __post_init__(self) -> None:
        """Validate parameters and derive the default anchor mismatch tolerance."""
        if self.max_mismatch >= self.min_anchor:
            raise ValueError(
                f"max_mismatch ({self.max_mismatch}) must be < min_anchor "
                f"({self.min_anchor}) to keep MUT-unique k-mers unambiguous"
            )
        if self.anchor_max_mismatch is None:
            object.__setattr__(self, "anchor_max_mismatch", max(self.max_mismatch, self.k // 12))


def decide_call(
    n_mut: int,
    vaf: float,
    n_mut_unique: int,
    evidence: str,
    min_mut_molecules: int,
    vaf_floor: float,
) -> str:
    """Return POSITIVE, NEGATIVE or UNRESOLVED from per-variant counts and thresholds."""
    if evidence == "kmer" and n_mut_unique == 0:
        return "UNRESOLVED"
    if n_mut >= min_mut_molecules and vaf >= vaf_floor:
        return "POSITIVE"
    return "NEGATIVE"


PRESETS = {
    "wgs": Params(input_level="wgs", max_mismatch=1),
    "min3": Params(input_level="min3", max_mismatch=1),
    "min1": Params(input_level="min1", max_mismatch=2),
}
