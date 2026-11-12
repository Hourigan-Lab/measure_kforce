# kforce — k-mer force-caller for MRD

Force-call known, patient-specific SV / ITDs / indels in deep sequencing data and decide their MRD (measurable-residual-disease) status. Instead of discovering variants, kforce asks, for each variant already in a manifest: is this exact junction still present in this sample, and at what fraction?

The core signal is alignment-free k-mers: an SV / ITD / indel junction creates sequence that exists in neither parent locus, so a k-mer spanning the junction is a specific signature of the mutant allele, robust to the soft-clipping / mismapping that defeats coordinate-based genotyping, and usable on UMI-deduplicated data.

---

## Install

```bash
git clone <repository-url> && cd measure_kforce
pixi install
pixi run test
```

Without pixi: Python ≥ 3.10 with `pysam`, `openpyxl`, `numpy`, `pandas` and (for the
tests) `pytest`; run from the repository root, e.g. `python -m kforce.cli call ...`.

## Inputs

| Input                                     | Notes                                                                                                             |
| ----------------------------------------- | ----------------------------------------------------------------------------------------------------------------- |
| **Alignment** (`--bam`)                   | Sorted, indexed BAM/CRAM. For UMI data, use the collapsed BAM (one record ≈ one molecule).                        |
| **Reference** (`--reference`)             | Indexed FASTA (`.fai`) of the genome the alignment was made against. A mismatch silently degrades classification. |
| **Manifest** (`--manifest`)               | `.xlsx` of the known variants, first worksheet, header row required (below).                                      |
| **Baseline** (`--baseline-bam`, optional) | A high-tumor (e.g. diagnosis) BAM/CRAM of the same patient; enables patient-specific references.                  |

Manifest columns (names configurable via the Python API, `load_manifest(columns=...)`):

```
Patient  Chromosome  Start  Stop  Alt  Ref  Category  VariantType  Group  T.VAF  N.VAF
```

`Chromosome` may be bare (`11`) or `chr`-prefixed contigs are normalized to `chr…`, so the
alignment and FASTA must use UCSC-style names. SVs are written as VCF-4.2 breakend ALT
strings (`G]chr9:20412833]`, `[chr11:118487577[AC`); indels use `Ref`/`Alt` with `-`
for the empty allele (a deletion has `Alt='-'`; an insertion has `Ref='-'` and `Stop=Start-1`);
`<DEL>`, `<INS>` and `<DUP:TANDEM>` are accepted as symbolic alleles. Only `Patient`,
`Chromosome`, `Start`, `Alt`, `Ref` are required; `T.VAF` is carried through for reference.

---

## How it works

Per variant (manifest row):

1. **Load** (`kforce/manifest.py`): read the manifest (xlsx), parse the variant (VCF-BND breakends, `Ref/Alt`/`-` indels,`<DUP:TANDEM>`, and the symbolic structural alts `<DEL>` (deletion of `[Start, Stop]`) and `<INS>` (insertion with no sequence given)), normalize contig names. SNVs are recognized and skipped (single-base ≈ sequencing error — the small-variant force-caller's job).
2. **Reconstruct** (`kforce/alleles.py`): build the local **WT** and **MUT** (junction)
   sequences and the **MUT-unique** k-mer set: k-mers in MUT absent from the reference and `> 2·D` mismatches from any local reference k-mer (so a within-`D` read match can only come from the mutant allele).
3. **Fetch** (`kforce/reads.py`): pull fragments over a window around each breakpoint (both loci for SVs), collapsing mates/splits to one molecule; a per-locus read cap bounds ultra-deep regions.
4. **Classify** (`kforce/classify.py`, `kforce/bothend.py`): label each fragment MUT / WT / ambiguous via complementary signals:
   - **MUT-unique k-mers** (exact or within `D` mismatches): the primary signal.
   - **Anchor-pair bridge**: catches junctions with untemplated bases / microhomology that shift the exact k-mers. For breakends, gated to far-apart loci (or different chromosomes) so a wild-type read can't co-occur with both anchors. For tandem duplications the geometry is safe at any span once the two k-mer anchors don't overlap (≥ 2×k): the bridge only fires when a read's reference position jumps backward from Stop to Start, which a wild-type read traversing the span in normal (forward) order can never show, so short-to-medium ITDs (e.g. FLT3-ITD) get this signal too, not just large PTDs.
   - **Both-end flanks**: for an UNRESOLVED junction in a repeat (no junction-unique k-mer), a fragment carrying unique k-mers from both flanks can only come from the rearranged molecule (alignment-free).
   - **Novel-insertion anchor-gap**: for symbolic `<INS>` (no sequence to build MUT k-mers from), a fragment carrying the two reference flank anchors separated by ≥ min_ins_gap untemplated bases is MUT; adjacent anchors (gap ≈ 0) are WT/reference. Works in both genome mode (anchor-gap) and baseline mode (observed insert with anchor-gap fallback).
   - **Discordant pairs** (opt-in, `--discordant-bothend`): for far-apart breakends that do have MUT-unique k-mers, also union the both-end signal so a read pair whose mates straddle the junction (no junction-spanning read) is counted MUT (see below).
5. **Call** (`kforce/call.py`): `VAF = n_mut/(n_mut+n_wt)`; **POSITIVE** if `n_mut ≥ min_mut_molecules AND vaf ≥ vaf_floor`. Statuses: `POSITIVE`, `NEGATIVE`,
   `UNRESOLVED` (no specific signal — e.g. a repeat junction with no reachable unique flank).
   The decision lives in one place (`config.decide_call`), shared by the live caller and `kforce recall` so they can't diverge. Breakend VAF is per-end: each breakpoint gives an independent estimate `n_mut/(n_mut+n_wt_end)` against its own reference population, and the reported `vaf` is their coverage-weighted average `2·n_mut/(2·n_mut+n_wt_a+n_wt_b)` (pooling both references against the single shared junction numerator would deflate it).
   The two per-breakpoint values are exposed as `vaf_bkpt_a` / `vaf_bkpt_b`. A large A-vs-B divergence flags an unreliable breakpoint (repeat / CNV / mapping artifact / one-sided).

### Patient-specific baseline reference (`--baseline-bam`)

The genome reference doesn't match a given patient (germline SNPs/indels in flanks; manifest indel-representation errors), costing MRD sensitivity. With a baseline (diagnosis,high-tumor) BAM, kforce rebuilds each variant's WT/MUT from the patient's own reads (`kforce/baseline.py`): germline-polished flanks (per-position consensus from the baseline pileup, `≥ --min-mapq`) + the junction/insert actually observed in the variant-supporting reads. For local indels, "observed" means comparing each baseline read's literal sequence over a fixed window around the declared edit (`observed_local_haplotype`, `--baseline-haplotype-margin`) against the manifest's own predicted sequence: trusting the manifest if enough reads match it exactly, otherwise trusting a single alternative only if it clearly dominates every competing one. For breakends/ITDs, "observed" means the anchor-pair bridge gap consensus (`observed_insert`). Calls are tagged `reference=genome` or
`reference=baseline`.

### Combined two-tier call (`--combined`, recommended for MRD)

Genome first, the (expensive) baseline only as a second-line check on the genome-negatives. Final call = genome if POSITIVE else baseline (per-variant `reference` column records which tier called it). Never loses a genome-positive while keeping every baseline recovery, and runs the baseline model only where needed.

### Tunable thresholds + re-thresholding without re-analysis

The positive call and all matching params are options on `kforce call`: `--min-mut-molecules` (N, default 2), `--vaf-floor` (default 0), `--min-mut-kmers` (3), `--min-wt-kmers` (3), `--max-mismatch` (D, 1), `--k` (25), `--window` (600), `--max-reads-per-locus` (0=unlimited), `--min-mapq` (20 — skip reads below this MAPQ everywhere reads are used as evidence), `--baseline-haplotype-margin` (20bp — `--baseline-bam` local-indel correction window, see `kforce/config.py` for the calibration notes), `--min-indel-size` (3); defaults in `kforce/config.py`. `--min-indel-size` filters the manifest at load: simple insertions/deletions whose size is ≤ N bp are not analyzed (default
3: only indels larger than 3 bp; 0 loads all). Size = the deleted span for a deletion, the inserted-base count for an insertion; symbolic `<INS>` (unknown size), structural `<DEL>` spans, and SV/breakend junctions are always kept. Because the call depends only on stored counts, the positive-call thresholds can be re-applied to existing outputs with no re-analysis via `kforce recall` / `kforce/recall.py` (the matching params `min-mut-kmers`/`max-mismatch`/ `k`/`min-mapq`/`baseline-haplotype-margin` change read fetching/classification and do need a rerun).

---

## Usage

```bash
# genome reference
python -m kforce.cli call --bam s.cram --manifest m.xlsx --reference hg38.fa --out out/

# patient baseline reference (germline-polished + observed junctions)
python -m kforce.cli call --bam mrd.cram --baseline-bam diagnosis.cram \
    --manifest m.xlsx --reference hg38.fa --emit-baseline-detail --out out/

# combined two-tier (recommended for MRD; loss-free)
python -m kforce.cli call --combined --bam mrd.cram --baseline-bam diagnosis.cram \
    --manifest m.xlsx --reference hg38.fa --out out/

# tighten the positive threshold
python -m kforce.cli call ... --min-mut-molecules 3 --vaf-floor 0.0001

# emit an IGV BAM of the supporting reads (KF tag), and recover discordant pairs
python -m kforce.cli call ... --emit-bam --discordant-bothend
```

Other subcommands:

- `kforce inspect --manifest m.xlsx --reference hg38.fa` : dump per-variant reconstructed MUT/WT sequences + MUT-unique k-mers for review (no BAM).
- `kforce recall --in <dir|tsv> --out <dir> --min-mut-molecules N [--vaf-floor F]` : re-derive calls from existing `*.kforce.tsv` with new thresholds, no re-analysis.
- `kforce aggregate '*/sample.kforce.tsv' --out cohort.tsv`: stack tidy TSVs.

`--min-indel-size N` skips small indels at load (see above). `--emit-reads` writes per-read evidence; `--emit-baseline-detail` writes the germline/indel corrections used.

**Per-sample outputs:** `<sample>.kforce.tsv` (tidy, one row/variant), `.summary.tsv` +
`.summary.html`, `.json`; `.reads.tsv` / `.baseline_detail.{tsv,txt}` / `.support.bam` when requested. Key tidy-TSV columns: `n_mut`, `n_wt`, `depth`, `vaf`, `vaf_bkpt_a` / `vaf_bkpt_b` (per-breakpoint VAF, breakends only), `call`, `reference` (`genome` | `baseline`, which tier made the call in `--combined` mode), `evidence` (`kmer` | `bothend` | `kmer+bothend` | `novel_ins`), `n_mut_unique`, `manifest_t_vaf`, `warnings`.

`--emit-bam` writes `<sample>.kforce.support.bam` (sorted + indexed): every MUT-supporting read, tagged `KF:Z:<category>` so IGV's *Color/Group by tag → KF* shows the evidence type — `KMER` (junction MUT-unique k-mers), `ANCHOR` (anchor-pair bridge), `BOTHEND` (both-end / discordant pair), `INS` (symbolic `<INS>` gap) — plus ~10% as many `WT` reads per variant for side-by-side comparison at the locus.

`--discordant-bothend` (opt-in, off by default): for far-apart breakends that already have MUT-unique k-mers, also count discordant read pairs (mates straddling the junction with no junction-spanning read) as MUT via the both-end signal. Requiring junction-spanning sequence is what makes the default robust to PE chimeric-ligation artifacts, so enabling this raises false-positive risk — consider raising `--min-mut-molecules` alongside it.

---

## Known limitations

- **Repeat-embedded junctions.** When both breakpoints sit in the same repeat family no junction-unique k-mer exists. Such variants are resolved through the both-end path where possible and otherwise reported `UNRESOLVED` (never a confident negative). The default `bothend_min_kmers=30` is deliberately conservative; lowering it trades specificity for recall at these loci.
- **Hypervariable tandem repeats.** Baseline correction can lock onto a consistent
  alternative repeat length that is a germline repeat polymorphism rather than the declared somatic variant. A baseline-vs-sample consistency check or a curated region exclusion list would be needed; neither is implemented.
- **Pre-existing structural artifacts** (segmental duplications, NUMTs) can recur across unrelated samples in the breakend / novel-insertion paths. A region blacklist is advisable for production use; none ships with kforce.
- **`--baseline-haplotype-margin`** is a fixed-width window: borderline loci near its edge can flip between corrected / uncorrected for nearby margin values.
- **Spanning-only breakends.** By default a breakend needs at least one junction-spanning read; a molecule evidenced only by a discordant pair is not counted. This keeps the caller robust to paired-end chimeric-ligation artifacts. `--discordant-bothend` recovers those molecules at higher false-positive risk (pair it with a higher `--min-mut-molecules`).
- **Scope.** SNVs are skipped; kforce is a force-caller for indels and structural variants.
- **Review low-molecule calls.** At MRD depth a couple of molecules can clear
  `min_mut_molecules`; inspect `--emit-bam` output in IGV before treating very-low-support positives as confident.

## Citation

If you use kforce, please cite the accompanying manuscript: *[citation to be added on
publication]*.

## License

See [`LICENSE`](LICENSE).
