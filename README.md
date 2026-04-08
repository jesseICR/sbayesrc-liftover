# SBayesRC snp.info hg19 to hg38 Liftover

Reproducible pipeline to lift the [SBayesRC](https://github.com/zhilizheng/SBayesRC) `snp.info` file (7,356,518 SNPs, autosomes 1-22) from GRCh37/hg19 to GRCh38/hg38.

## Quick Start

```bash
# Requirements: Python 3.10+, curl
bash main.sh
```

That's it. The script creates a virtual environment in `tools/venv/`, installs dependencies from `requirements.txt`, downloads `snp.info` and all reference files on first run (see [Storage](#storage)), caches everything, and produces the output files.

**Nothing is installed globally.** All Python packages, binaries, and reference data are stored inside the repo directory (`tools/`, `tmp/`). Delete the repo and there is zero trace left on your system.

Pre-built output files are also available as [GitHub Release](https://github.com/jesseICR/sbayesrc-liftover/releases) assets -- no pipeline run required.

### Docker

Pull the pre-built image from GHCR:

```bash
docker pull ghcr.io/jesseicr/sbayesrc-liftover:latest
```

Or build locally:

```bash
docker build -t sbayesrc-liftover .
```

Run (mount the repo directory so downloads persist and output is accessible):

```bash
docker run --rm -v $(pwd):/data ghcr.io/jesseicr/sbayesrc-liftover:latest
```

The ~5 GB of downloaded reference files are cached in `tools/` and reused on subsequent runs.

## Why This Exists

UCSC liftOver works well for the vast majority of SNPs, but it makes errors in **segmental duplication regions** (NBPF genes on chr1, HLA on chr6, etc.) where it maps variants to paralogous copies at the wrong position. These errors can be invisible to simple allele-vs-reference checks when the allele coincidentally matches the reference at the wrong location.

Rather than attempt to override these errors (which introduces its own risks), this pipeline takes a **conservative approach**: every SNP in the final output must be independently confirmed by both UCSC liftOver (or dbSNP rescue) and dbSNP, and validated against the hg38 reference FASTA. Any disagreement results in exclusion.

## Inclusion Criteria

Every SNP in the final output (`sbayesrc_hg38.csv`) must satisfy **all three**:

1. **dbSNP match** -- the rsID exists in the current dbSNP release with a matching hg38 chromosome and position
2. **liftOver agreement** -- if UCSC liftOver succeeded, its hg38 position must equal the dbSNP position. (If liftOver failed but dbSNP has a position, the SNP is included as a "rescue".)
3. **FASTA confirmation** -- the reference allele that dbSNP reports for the rsID must match the actual hg38 FASTA reference base at that position, and at least one of the SNP's alleles (A1 or A2, accounting for strand) must also match the reference

## Pipeline

### Step 1: UCSC liftOver

Converts hg19 coordinates to hg38 using the UCSC liftOver tool and the hg19-to-hg38 chain file. Writes a 6-column BED file to detect strand flips from the chain alignment. SNPs that map to a different chromosome are treated as unmapped.

**Output:** `pos_liftover` for each SNP -- the hg38 position that liftOver produced, or -1 if liftOver failed.

### Step 2: dbSNP cross-validation

Compares liftOver positions against dbSNP (GRCh38) using rsID as the join key. Every SNP is classified into one of five statuses:

| Status | Condition | Included in final output? |
|--------|-----------|:-------------------------:|
| `confirmed` | liftOver succeeded AND dbSNP agrees on chrom + position | Yes |
| `rescue` | liftOver failed AND dbSNP has a position for this rsID | Yes (pending FASTA check) |
| `conflict` | liftOver succeeded AND dbSNP has a **different** position | **No** |
| `no_dbsnp` | rsID not found in the current dbSNP release | **No** |
| `unmapped` | liftOver failed AND no dbSNP entry | **No** |

Conflicts are the key safety mechanism: if UCSC and dbSNP disagree on where a SNP is, neither is trusted.

### Step 3: FASTA reference validation

For all `confirmed` and `rescue` SNPs, fetches the hg38 FASTA reference base at the hg38 position and checks:

1. **dbSNP ref vs FASTA ref** -- does the reference allele that dbSNP reports match the actual FASTA base? If not, the SNP is reclassified as `fasta_mismatch` and excluded.
2. **Allele vs FASTA ref** -- does at least one of the SNP's alleles (A1 or A2, or their strand complements) match the FASTA reference? If not, the SNP is reclassified as `allele_mismatch` and excluded.

### Step 4: Allele annotation

For all SNPs that passed steps 1-3:

- **Strand-flipped liftOver SNPs:** alleles are complemented (A<->T, C<->G) based on the strand flip detected from the chain alignment in step 1
- **Rescue SNPs:** alleles are complemented if needed to match the FASTA reference strand
- **Ref/alt determination:** determines whether A1 or A2 matches the hg38 reference

## Results

*Numbers below are from the most recent pipeline run. Re-run `bash main.sh` to regenerate.*

| Status | Count | In final output? |
|--------|------:|:-----------------:|
| `confirmed` | 7,352,740 | Yes |
| `rescue` | 2,213 | Yes |
| `conflict` | 25 | No |
| `no_dbsnp` | 1,538 | No |
| `unmapped` | 2 | No |
| `fasta_mismatch` | 0 | No |
| `allele_mismatch` | 0 | No |
| **Total input** | **7,356,518** | |
| **Total in sbayesrc_hg38.csv** | **7,354,953** | |

- **25 conflicts** between UCSC liftOver and dbSNP, all in segmental duplication and HLA regions (chr1 NBPF, chr3, chr6 HLA, chr11). These are discarded.
- **2,213 rescues** where liftOver failed but dbSNP provides a confirmed position.
- **1,538 SNPs** have no dbSNP entry (0.02%) and are excluded because they cannot be independently confirmed.
- **2 unmapped SNPs**: rs117553620 (chr17) and rs140636911 (chr19) -- no source resolves them.
- **0 FASTA mismatches** -- every included SNP has dbSNP ref matching the hg38 FASTA reference base.
- **12 strand-flipped SNPs** complemented (9 from liftOver chain alignment, 3 rescue SNPs).

## Logging

Every pipeline run writes a timestamped log file to `logs/` (e.g., `logs/run_20260408_161500.log`). The log captures all console output including per-step SNP counts and the final summary. The `logs/` directory is gitignored.

## Input

The canonical SBayesRC `snp.info` file (tab-delimited):

| Column | Description |
|--------|-------------|
| Chrom | Chromosome (1-22) |
| ID | rsID |
| Index | SBayesRC internal index |
| GenPos | Genetic position |
| PhysPos | Physical position (hg19) |
| A1 | Effect allele |
| A2 | Other allele |
| A1Freq | Effect allele frequency |
| N | Sample size |
| Block | LD block |

## Output

Two output files are produced:

### `sbayesrc_hg38.csv` (primary output)

Clean, minimal file with one row per SNP that passed all three inclusion criteria:

| Column | Description |
|--------|-------------|
| chrom | Chromosome |
| pos | hg38 position |
| ref | Reference allele (hg38) |
| alt | Alternate allele (hg38) |
| rsid | rsID |

### `sbayesrc_liftover_results.csv` (verbose results)

Full liftover details for all 7,356,518 SNPs (including excluded):

| Column | Description |
|--------|-------------|
| chrom | Chromosome |
| ID | rsID |
| pos_hg19 | Original hg19 position |
| pos_hg38 | Final hg38 position (-1 if excluded) |
| pos_liftover | Position from UCSC liftOver (-1 if liftOver failed) |
| pos_dbsnp | Position from dbSNP (-1 if rsID not in dbSNP) |
| A1, A2 | Original alleles from snp.info |
| A1_hg38, A2_hg38 | Alleles on the hg38 strand (complemented if strand-flipped) |
| dbsnp_ref | Reference allele reported by dbSNP (empty if no dbSNP entry) |
| fasta_ref | Actual hg38 FASTA reference base at pos_hg38 (empty if excluded) |
| ref_match | Which allele matches ref: `A1_hg38` or `A2_hg38` (empty if excluded) |
| strand_flip | Whether alleles were complemented |
| status | SNP classification (see [Step 2](#step-2-dbsnp-cross-validation) and [Step 3](#step-3-fasta-reference-validation)) |
| Index, GenPos, A1Freq, N, Block | Preserved from snp.info |

## Runtime and Storage

**Benchmarked on:** 128 cores, 503 GiB RAM (AMD Threadripper PRO 5995WX)

### Per-Step Runtime

| Step | Description | Runtime |
|------|-------------|---------|
| 0 | Setup (download tools, FASTA, chain file, dbSNP VCF) | ~30 min first run, <1 s cached |
| 0b | Build dbSNP lookup from VCF (first run only) | ~15-30 min |
| 1 | UCSC liftOver (7.36M SNPs) | ~40 s |
| 2 | dbSNP cross-validation | ~10 s |
| 3 | FASTA reference validation | ~10 s |
| 4 | Allele annotation | ~10 s |
| -- | Write output CSVs | ~20 s |

**Total pipeline runtime: ~100 seconds** (cached run, all reference files present).

First-run download time depends on network speed. The dbSNP VCF (~28 GB) and hg38 FASTA (~900 MB) are the largest downloads. The dbSNP lookup extraction streams the full VCF once (~15-30 min) and is cached for subsequent runs.

### Storage

| Directory / File | Size | Contents |
|------------------|-----:|---------|
| `tools/hg38.fa` | 3.1 GB | GRCh38 reference FASTA (UCSC, decompressed) |
| `tools/hg38.fa.gz` | 939 MB | Compressed FASTA (kept for re-extraction) |
| `tools/dbsnp_lookup.tsv` | 185 MB | rsID-to-position table (streamed from 28 GB dbSNP VCF) |
| `tools/venv/` | ~100 MB | Python virtual environment |
| `tools/bin/liftOver` | 24 MB | UCSC liftOver binary |
| `tools/hg19ToHg38.over.chain.gz` | 224 KB | hg19-to-hg38 chain file |
| `tmp/` | ~250 MB | Intermediate BED files |
| **Total** | **~5 GB** | |

On first run, the dbSNP VCF (~28 GB) is **streamed directly** from NCBI FTP and decompressed on the fly -- only the small lookup TSV (185 MB) is saved to disk. The hg38 FASTA (~900 MB compressed) is downloaded and decompressed. All files are cached in `tools/` and reused on subsequent runs.

## Directory Structure

```
.
├── main.sh                     Entry point (venv setup + logging + run)
├── liftover.py                 Pipeline (single file, self-contained)
├── requirements.txt            Python dependencies (pandas, pysam)
├── snp.info                    Input: SBayesRC SNPs, hg19 (downloaded from GitHub Release)
├── Dockerfile
├── .dockerignore
├── .gitignore
├── .github/
│   └── workflows/
│       └── docker-publish.yml  Builds and publishes Docker image to GHCR
├── README.md
├── logs/                       Timestamped run logs (gitignored)
├── tools/                      Downloaded reference files (gitignored)
│   ├── bin/liftOver            UCSC liftOver binary
│   ├── venv/                   Python virtual environment
│   ├── hg19ToHg38.over.chain.gz
│   ├── hg38.fa                 GRCh38 reference FASTA
│   ├── hg38.fa.fai             FASTA index (created by pysam)
│   └── dbsnp_lookup.tsv        rsID-to-position table (streamed from dbSNP VCF)
└── tmp/                        Intermediate files (gitignored)
    └── *.bed                   Temporary liftOver files
```

## Requirements

- Python 3.10+ (a virtual environment is created automatically)
- `curl` (for dbSNP VCF download)

## Known Limitations

- **~1,540 SNPs have no dbSNP entry** in the current release (0.02% of 7.36M). These are excluded from the final output (`no_dbsnp` status) because they cannot be independently confirmed.
- **Strand-ambiguous SNPs (A/T, C/G):** Strand flips are detected from the chain alignment, not from alleles. For A/T and C/G SNPs, complementing swaps labels but not nucleotides. This is correct but may look surprising in the output.
- **dbSNP `latest_release`:** The pipeline streams from NCBI's `latest_release` URL, so the dbSNP version changes when NCBI publishes updates. The cached `dbsnp_lookup.tsv` is not automatically refreshed -- delete it to re-stream from the latest release.

## Data Sources

- **UCSC liftOver**: https://hgdownload.soe.ucsc.edu/admin/exe/
- **hg19-to-hg38 chain file**: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
- **GRCh38 reference FASTA**: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
- **dbSNP VCF**: https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/
- **SBayesRC**: https://github.com/zhilizheng/SBayesRC
