"""Write a small IGV-ready BAM of the reads supporting each MUT call."""
from __future__ import annotations
import os
import random
import pysam
from .config import Params


def write_support_bam(
    bam: pysam.AlignmentFile,
    support: list[dict],
    out_path: str,
    params: Params,
    wt_fraction: float = 0.10,
    seed: int = 0,
) -> dict:
    """Write and index a BAM of tagged MUT and sampled WT reads."""
    valid = set(bam.references)
    W = params.window
    tmp = out_path + ".unsorted.bam"
    written: set[tuple] = set()
    n_mut_reads = n_wt_reads = 0

    with pysam.AlignmentFile(tmp, "wb", header=bam.header) as out:
        for entry in support:
            qtag = dict(entry["mut"])
            n_keep = round(wt_fraction * len(entry["mut"]))
            wt = sorted(entry["wt"])
            if n_keep and wt:
                keep = wt if n_keep >= len(wt) else random.Random(seed).sample(wt, n_keep)
                for qn in keep:
                    qtag.setdefault(qn, "WT")
            if not qtag:
                continue
            for chrom, pos in entry["loci"]:
                if chrom not in valid:
                    continue
                try:
                    reads = list(bam.fetch(chrom, max(0, pos - 1 - W), pos - 1 + W))
                except (OSError, ValueError):
                    continue
                for r in reads:
                    if r.is_secondary:
                        continue
                    cat = qtag.get(r.query_name)
                    if cat is None:
                        continue
                    key = (r.query_name, r.flag, r.reference_id, r.reference_start)
                    if key in written:
                        continue
                    written.add(key)
                    r.set_tag("KF", cat, value_type="Z")
                    out.write(r)
                    if cat == "WT":
                        n_wt_reads += 1
                    else:
                        n_mut_reads += 1

    pysam.sort("-o", out_path, tmp)
    pysam.index(out_path)
    os.remove(tmp)
    return {
        "bam": out_path,
        "n_mut_reads": n_mut_reads,
        "n_wt_reads": n_wt_reads,
        "n_variants": sum(1 for e in support if e["mut"]),
    }
