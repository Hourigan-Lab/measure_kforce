"""Small dependency-free sequence helpers."""
from __future__ import annotations
from typing import Iterator

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    """Return the reverse complement of a sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


def iter_kmers(seq: str, k: int) -> Iterator[str]:
    """Yield every length-k substring of a sequence."""
    if k <= 0:
        raise ValueError("k must be positive")
    for i in range(len(seq) - k + 1):
        yield seq[i : i + k]


def hamming(a: str, b: str) -> int:
    """Return the Hamming distance between two equal-length strings."""
    if len(a) != len(b):
        raise ValueError(f"hamming() needs equal lengths, got {len(a)} and {len(b)}")
    return sum(1 for x, y in zip(a, b) if x != y)


def min_hamming(query: str, targets) -> int:
    """Return the minimum Hamming distance from a query to any target."""
    best = len(query) + 1
    for t in targets:
        d = 0
        for x, y in zip(query, t):
            if x != y:
                d += 1
                if d >= best:
                    break
        if d < best:
            best = d
            if best == 0:
                break
    return best


def is_acgt(seq: str) -> bool:
    """Return True if a sequence contains only A, C, G and T."""
    return bool(seq) and set(seq) <= {"A", "C", "G", "T"}
