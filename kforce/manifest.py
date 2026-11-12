"""Load the variant manifest into normalized Variant records."""
from __future__ import annotations
import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional
import openpyxl


class BndOrient(str, Enum):
    """Orientation of a breakend (which side is kept and whether the mate is flipped)."""
    LEFT_FWD = "t[p["
    LEFT_REV = "t]p]"
    RIGHT_FWD = "]p]t"
    RIGHT_REV = "[p[t"


@dataclass(frozen=True)
class Breakend:
    """Parsed VCF breakend ALT: mate locus, orientation, reference base and insert."""
    mate_chrom: str
    mate_pos: int
    orient: BndOrient
    ref_base: str
    insert: str


@dataclass(frozen=True)
class Variant:
    """One manifest row as a normalized variant."""
    patient: str
    chrom: str
    pos: int
    stop: int
    ref: str
    alt: str
    category: str
    variant_type: str
    group: str
    t_vaf: Optional[float]
    n_vaf: Optional[float]
    breakend: Optional[Breakend]

    @property
    def kind(self) -> str:
        """Return the coarse variant class that selects the reconstruction strategy."""
        if self.breakend is not None:
            return "breakend"
        alt = self.alt.upper()
        if alt.startswith("<DUP"):
            return "tandem_dup"
        if alt.startswith("<DEL"):
            return "deletion"
        if alt.startswith("<INS"):
            return "insertion"
        if alt.startswith("<") or alt.startswith("."):
            return "other"
        ref = "" if self.ref in ("", "-") else self.ref.upper()
        a = "" if alt in ("", "-") else alt
        if ref and a and len(ref) == 1 and len(a) == 1:
            return "snv"
        if not ref and a:
            return "insertion"
        if ref and not a:
            return "deletion"
        if ref and a:
            return "mnv" if len(ref) == len(a) else "delins"
        return "other"

    @property
    def indel_size(self) -> Optional[int]:
        """Return the size of a simple indel in bp, or None if not a sized indel."""
        k = self.kind
        if k == "deletion":
            return self.stop - self.pos + 1
        if k == "insertion":
            a = self.alt.upper()
            if not a or a == "-" or a.startswith("<"):
                return None
            return len(a)
        return None

    @property
    def variant_id(self) -> str:
        """Return a stable key of the form chrom:pos:type:mate."""
        vt = self.variant_type.replace(" ", "_")
        if self.breakend is not None:
            mate = f"{self.breakend.mate_chrom}:{self.breakend.mate_pos}"
        else:
            mate = str(self.stop)
        return f"{self.chrom}:{self.pos}:{vt}:{mate}"


def normalize_contig(name: str) -> str:
    """Return a chr-prefixed contig name."""
    name = str(name).strip()
    if name.lower().startswith("chr"):
        return "chr" + name[3:]
    return "chr" + name


_RE_LEFT = re.compile(r"^([ACGTNacgtn]+)([\[\]])([^:\[\]]+):(\d+)\2$")
_RE_RIGHT = re.compile(r"^([\[\]])([^:\[\]]+):(\d+)\1([ACGTNacgtn]+)$")


def parse_breakend(alt: str) -> Breakend:
    """Parse a VCF-4.2 breakend ALT string into a Breakend."""
    alt = alt.strip()

    m = _RE_LEFT.match(alt)
    if m:
        t, bracket, mate_chrom, mate_pos = m.groups()
        orient = BndOrient.LEFT_FWD if bracket == "[" else BndOrient.LEFT_REV
        ref_base = t[0].upper()
        insert = t[1:].upper()
        return Breakend(normalize_contig(mate_chrom), int(mate_pos), orient, ref_base, insert)

    m = _RE_RIGHT.match(alt)
    if m:
        bracket, mate_chrom, mate_pos, t = m.groups()
        orient = BndOrient.RIGHT_FWD if bracket == "]" else BndOrient.RIGHT_REV
        ref_base = t[-1].upper()
        insert = t[:-1].upper()
        return Breakend(normalize_contig(mate_chrom), int(mate_pos), orient, ref_base, insert)

    raise ValueError(f"not a breakend ALT: {alt!r}")


def is_breakend(alt: str) -> bool:
    """Return True if an ALT string is in breakend notation."""
    return bool(_RE_LEFT.match(alt.strip()) or _RE_RIGHT.match(alt.strip()))


DEFAULT_COLUMNS = {
    "patient": "Patient",
    "chrom": "Chromosome",
    "pos": "Start",
    "stop": "Stop",
    "ref": "Ref",
    "alt": "Alt",
    "category": "Category",
    "variant_type": "VariantType",
    "group": "Group",
    "t_vaf": "T.VAF",
    "n_vaf": "N.VAF",
}


def _to_float(v) -> Optional[float]:
    """Convert a cell value to float, returning None if empty or invalid."""
    if v is None or v == "":
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def filter_by_indel_size(variants: list[Variant], min_size: int) -> tuple[list[Variant], int]:
    """Drop simple indels of at most a given size and return (kept, number dropped)."""
    if min_size <= 0:
        return list(variants), 0
    kept, dropped = [], 0
    for v in variants:
        s = v.indel_size
        if s is not None and s <= min_size:
            dropped += 1
        else:
            kept.append(v)
    return kept, dropped


def load_manifest(
    path, patient: Optional[str] = None, columns: Optional[dict] = None, sheet: Optional[str] = None
) -> list[Variant]:
    """Read an .xlsx manifest into a list of Variant records."""
    columns = {**DEFAULT_COLUMNS, **(columns or {})}
    wb = openpyxl.load_workbook(Path(path), read_only=True, data_only=True)
    ws = wb[sheet] if sheet else wb.worksheets[0]
    if hasattr(ws, "reset_dimensions"):
        ws.reset_dimensions()

    rows = list(ws.iter_rows(values_only=True))
    if not rows:
        return []
    header = [str(h).strip() if h is not None else "" for h in rows[0]]
    idx = {field: header.index(col) for field, col in columns.items() if col in header}

    required = ("chrom", "pos", "alt", "ref")
    missing = [columns[f] for f in required if f not in idx]
    if missing:
        raise ValueError(f"manifest {path} missing required columns: {missing}")

    def get(row, field):
        """Return a row's value for a named field, or None if absent."""
        i = idx.get(field)
        return row[i] if i is not None and i < len(row) else None

    variants: list[Variant] = []
    for row in rows[1:]:
        if row is None or all(c is None for c in row):
            continue
        pat = get(row, "patient")
        if patient is not None and pat is not None and str(pat) != patient:
            continue
        alt = str(get(row, "alt") or "").strip()
        bnd = parse_breakend(alt) if is_breakend(alt) else None
        pos = int(get(row, "pos"))
        stop_raw = get(row, "stop")
        stop = int(stop_raw) if stop_raw not in (None, "") else pos
        variants.append(
            Variant(
                patient=str(pat) if pat is not None else "",
                chrom=normalize_contig(get(row, "chrom")),
                pos=pos,
                stop=stop,
                ref=str(get(row, "ref") or "").strip().upper(),
                alt=alt,
                category=str(get(row, "category") or "").strip(),
                variant_type=str(get(row, "variant_type") or "").strip(),
                group=str(get(row, "group") or "").strip(),
                t_vaf=_to_float(get(row, "t_vaf")),
                n_vaf=_to_float(get(row, "n_vaf")),
                breakend=bnd,
            )
        )
    return variants
