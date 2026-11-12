"""Fetch sequenced fragments (molecules) overlapping a variant's loci."""
from __future__ import annotations
from typing import Iterable, Optional
import pysam
from .config import Params
from .manifest import Variant


def open_alignment(
    path: str, reference: Optional[str] = None, threads: int = 1
) -> pysam.AlignmentFile:
    """Open a BAM or CRAM file."""
    if str(path).endswith(".cram"):
        return pysam.AlignmentFile(path, "rc", reference_filename=reference, threads=threads)
    return pysam.AlignmentFile(path, "rb", threads=threads)


def variant_loci(v: Variant) -> list[tuple[str, int]]:
    """Return the genomic loci to fetch reads from for a variant."""
    loci = [(v.chrom, v.pos)]
    if v.breakend is not None:
        loci.append((v.breakend.mate_chrom, v.breakend.mate_pos))
    elif v.stop != v.pos:
        loci.append((v.chrom, v.stop))
    return loci


def fetch_fragments(
    bam: pysam.AlignmentFile,
    loci: Iterable[tuple[str, int]],
    params: Params,
) -> dict[str, list[str]]:
    """Return read sequences per fragment name over the given loci."""
    W = params.window
    cap = params.max_reads_per_locus
    by_end: dict[str, dict[bool, str]] = {}
    valid = set(bam.references)
    for chrom, pos in loci:
        if chrom not in valid:
            continue
        start = max(0, pos - 1 - W)
        end = pos - 1 + W
        scanned = 0
        try:
            for r in bam.fetch(chrom, start, end):
                if cap and scanned >= cap:
                    break
                scanned += 1
                if r.is_secondary or r.is_supplementary:
                    continue
                if r.mapping_quality < params.min_mapq:
                    continue
                seq = r.query_sequence
                if not seq:
                    continue
                by_end.setdefault(r.query_name, {})[bool(r.is_read1)] = seq.upper()
        except (OSError, ValueError):
            continue
    return {qn: list(ends.values()) for qn, ends in by_end.items()}
