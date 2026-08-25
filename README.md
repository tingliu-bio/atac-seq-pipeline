# ATAC-seq Analysis Workflow

This repo documents an ATAC-seq chromatin accessibility analysis pipeline. Core alignment and peak calling steps follow the [Lupien Lab chromatin accessibility workflow](https://github.com/LupienLab/pipeline-chromatin-accessibility), with modifications described below.

---

## 1. Alignment and Peak Calling

Core steps follow an established ATAC-seq pipeline (read QC/trimming, Bowtie2 alignment, duplicate marking, and MACS2 peak calling), configured via YAML.

**Custom reference genome:** Alignment was performed against a merged reference genome (hg38 + additional sequences) rather than standard hg38, to account for project-specific requirements.

The full alignment-to-peaks flow:
FASTQ → QC/trimming → Bowtie2 (hg38+HTLV1) → deduplication → MACS2 → filtered.narrowPeak

> **Note for IDR:** IDR uses `peaks.filtered.narrowPeak` generated from `filtered.dedup.sorted.bam` as input. If IDR is skipped, `peaks.filtered.merged.narrowPeak` may be used instead.

---

## 2. IDR Analysis

IDR was run per condition using two biological replicates to generate a reproducible peak set before cross-condition merging.

<pre>
Replicate 1 ──┐
              |→ IDR → condition-level reproducible peaks
Replicate 2 ──┘
</pre>

**Input files:** Use `*_peaks.filtered.narrowPeak` — these preserve the full 10-column narrowPeak format that IDR expects.

> **Do NOT use** `*_peaks.filtered.merged.narrowPeak`. These have been through `bedtools merge`, which collapses the column structure down to ~7 columns and breaks IDR scoring.

**IDR threshold:** We kept peaks at IDR ≤ 0.05, set via `--soft-idr-threshold 0.05`.

---

## 3. Consensus Peak Universe

Per-condition IDR peaks were merged into a single consensus universe for downstream analysis (DiffBind, DESeq2, annotation).

<pre>
 Condition 1 IDR peaks ─┐
 Condition 2 IDR peaks ─┼→ consensus_IDR_peaks.bed
 Condition 3 IDR peaks ─┘
</pre>

`consensus_IDR_peaks.bed` captures all reproducible open chromatin regions across conditions and serves as the shared coordinate space for downstream work.

Depending on your analysis, you may want condition-specific IDR peaks instead of (or in addition to) the consensus set.