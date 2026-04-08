# SBayesRC snp.info hg19 to hg38 Liftover

Reproducible pipeline to lift the [SBayesRC](https://github.com/zhilizheng/SBayesRC) `snp.info` file (7,356,518 SNPs, autosomes 1-22) from GRCh37/hg19 to GRCh38/hg38.

## Quick Start

```bash
# Requirements: Python 3.10+, curl
bash main.sh
```

That's it. The script creates a virtual environment in `tools/venv/`, installs dependencies from `requirements.txt`, downloads `snp.info` and all reference files on first run (see [Storage](#storage)), caches everything, and produces the output files.

**Nothing is installed globally.** All Python packages, binaries, and reference data are stored inside the repo directory (`tools/`, `tmp/`). Delete the repo and there is zero trace left on your system.

Pre-built output files are also available as [GitHub Release](https://github.com/jesseICR/sbayesrc-liftover/releases) assets — no pipeline run required.

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

UCSC liftOver works well for the vast majority of SNPs, but it makes errors in **segmental duplication regions** (NBPF genes on chr1, HLA on chr6, etc.) where it maps variants to paralogous copies at the wrong position. These errors are invisible to simple allele-vs-reference checks when the allele coincidentally matches the reference at the wrong location.

This pipeline cross-validates every liftOver position against dbSNP and the Ensembl REST API to catch and correct these errors. The result is a high-confidence hg38 coordinate set with ref/alt allele annotation for all 7.36M SBayesRC SNPs.

## Pipeline

| Step | Source | What it does |
|------|--------|--------------|
| 1 | UCSC liftOver | Position-based coordinate conversion using hg19-to-hg38 chain file. Writes 6-column BED to detect strand flips. |
| 2 | hg38 FASTA | Checks the reference allele at each liftOver position. Flags "suspects" where neither A1 nor A2 matches (possible wrong mapping). |
| 3 | dbSNP VCF | Cross-validates liftOver positions against [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) (GRCh38). Overrides disagreements, fills liftOver failures. |
| 4 | Ensembl API | Queries the [Ensembl REST API](https://rest.ensembl.org/) for remaining unmapped SNPs, unresolved suspects, and any unverified liftOver positions in the cache. Overrides when Ensembl disagrees with liftOver. |
| 5 | hg38 FASTA | Final allele annotation. Complements strand-flipped alleles, determines which allele matches the hg38 reference. |

### Cascade Logic

For each SNP, the position source is chosen by priority:

1. If **dbSNP disagrees** with liftOver, use dbSNP (rsID-based mapping is authoritative for segdup regions).
2. If **Ensembl disagrees** with liftOver (and liftOver wasn't already validated by dbSNP), use Ensembl.
3. If both dbSNP and Ensembl are unavailable for a SNP, **trust liftOver** (correct >99.99% of the time).
4. If liftOver failed, use **dbSNP or Ensembl as fallback**.
5. If no source can map the SNP, mark as **unmapped**.

## Results

| Category | Count |
|----------|------:|
| Mapped via liftOver (confirmed by dbSNP) | 7,354,279 |
| Mapped via dbSNP (override or rescue) | 2,237 |
| Unmapped (no source resolves) | 2 |
| **Total** | **7,356,518** |

- **24 liftOver errors corrected** by dbSNP in segmental duplication and HLA regions (1 additional disagreement skipped due to ref mismatch at the dbSNP position)
- **20 strand-flipped SNPs** complemented (14 detected by liftOver, 6 by FASTA ref check)
- **0 alleles with "neither" ref match** -- every mapped SNP has A1 or A2 matching the hg38 reference
- **2 unmapped SNPs**: rs117553620 (chr17) and rs140636911 (chr19)
- **dbSNP coverage: 99.98%** (7,354,979 / 7,356,518 rsIDs found in the current dbSNP release)

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

Clean, minimal file with one row per successfully mapped SNP:

| Column | Description |
|--------|-------------|
| chrom | Chromosome |
| pos | hg38 position |
| ref | Reference allele (hg38) |
| alt | Alternate allele (hg38) |
| rsid | rsID |

### `sbayesrc_liftover_results.csv` (verbose results)

Full liftover details for all 7,356,518 SNPs (including unmapped):

| Column | Description |
|--------|-------------|
| chrom | Chromosome |
| ID | rsID |
| pos_hg19 | Original hg19 position |
| pos_hg38 | Lifted hg38 position (-1 if unmapped) |
| A1, A2 | Original alleles from snp.info |
| A1_hg38, A2_hg38 | Alleles on the hg38 strand (complemented if strand-flipped) |
| ref_hg38 | Reference allele at the hg38 position |
| ref_match | Which allele matches ref: `A1_hg38`, `A2_hg38`, or `neither` |
| strand_flip | Whether alleles were complemented (strand flip detected) |
| method | Position source: `liftover`, `dbsnp`, `ensembl`, or empty (unmapped) |
| Index, GenPos, A1Freq, N, Block | Preserved from snp.info |

Both files are also checked into `outputs/` so results are available directly from the repository without running the pipeline.

## Exhaustive Ensembl Validation

By default, the pipeline queries Ensembl only for high-priority SNPs (unmapped + FASTA suspects). To cross-validate all liftOver positions that lack dbSNP confirmation:

```bash
ENSEMBL_FULL=1 bash main.sh
```

This queries remaining rsIDs through the Ensembl REST API (~3-4 hours, cached in `tmp/ensembl_cache.json`). Subsequent runs apply cached results instantly.

## Runtime and Storage

**Benchmarked on:** 128 cores, 503 GiB RAM (AMD Threadripper PRO 5995WX)

### Per-Step Runtime

| Step | Description | Runtime |
|------|-------------|---------|
| 0 | Setup (download tools, FASTA, chain file, dbSNP VCF) | ~30 min first run, <1 s cached |
| 0b | Build dbSNP lookup from VCF (first run only) | ~15-30 min |
| 1 | UCSC liftOver (7.36M SNPs) | 37 s |
| 2 | FASTA ref check (flag suspects) | 9 s |
| 3 | dbSNP cross-validation | 9 s |
| 4 | Ensembl cross-validation (high-priority only) | 12 s |
| 5 | Allele annotation (FASTA query + complement) | 7 s |
| -- | Write output CSVs | 19 s |

**Total pipeline runtime: ~100 seconds** (cached run, all reference files present).

First-run download time depends on network speed. The dbSNP VCF (~28 GB) and hg38 FASTA (~900 MB) are the largest downloads. The dbSNP lookup extraction streams the full VCF once (~15-30 min) and is cached for subsequent runs.

With `ENSEMBL_FULL=1`, Step 4 queries additional rsIDs via the Ensembl REST API, adding ~3-4 hours (also cached).

### Storage

| Directory / File | Size | Contents |
|------------------|-----:|---------|
| `tools/hg38.fa` | 3.1 GB | GRCh38 reference FASTA (UCSC, decompressed) |
| `tools/hg38.fa.gz` | 939 MB | Compressed FASTA (kept for re-extraction) |
| `tools/dbsnp_lookup.tsv` | 185 MB | rsID-to-position table (streamed from 28 GB dbSNP VCF) |
| `tools/venv/` | ~100 MB | Python virtual environment |
| `tools/bin/liftOver` | 24 MB | UCSC liftOver binary |
| `tools/hg19ToHg38.over.chain.gz` | 224 KB | hg19-to-hg38 chain file |
| `tmp/` | 502 MB | Intermediate BED files + Ensembl cache |
| **Total** | **~5 GB** | |

On first run, the dbSNP VCF (~28 GB) is **streamed directly** from NCBI FTP and decompressed on the fly — only the small lookup TSV (185 MB) is saved to disk. The hg38 FASTA (~900 MB compressed) is downloaded and decompressed. All files are cached in `tools/` and reused on subsequent runs.

## Directory Structure

```
.
├── main.sh                     Entry point (venv setup + run)
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
├── outputs/                    Pipeline outputs (available as GitHub Release assets)
│   ├── sbayesrc_hg38.csv       Clean hg38 coordinates (chrom, pos, ref, alt, rsid)
│   └── sbayesrc_liftover_results.csv  Verbose liftover details
├── tools/                      Downloaded reference files (gitignored)
│   ├── bin/liftOver            UCSC liftOver binary
│   ├── venv/                   Python virtual environment
│   ├── hg19ToHg38.over.chain.gz
│   ├── hg38.fa                 GRCh38 reference FASTA
│   ├── hg38.fa.fai             FASTA index (created by pysam)
│   └── dbsnp_lookup.tsv        rsID-to-position table (streamed from dbSNP VCF)
└── tmp/                        Intermediate files (gitignored)
    ├── ensembl_cache.json      Cached Ensembl API results
    └── *.bed, *.txt            Temporary liftOver files
```

## Requirements

- Python 3.10+ (a virtual environment is created automatically)
- `curl` (for dbSNP VCF download)

## Known Limitations

- **~1,540 SNPs have no dbSNP entry** in the current release (0.02% of 7.36M). These rely solely on liftOver positions. Use `ENSEMBL_FULL=1` to cross-validate them via the Ensembl API.
- **Strand-ambiguous SNPs (A/T, C/G):** Strand flips are detected from the chain alignment, not from alleles. For A/T and C/G SNPs, complementing swaps labels but not nucleotides. This is correct but may look surprising in the output.
- **2 unmapped SNPs** cannot be resolved by any source (liftOver, dbSNP, or Ensembl).

## Data Sources

- **UCSC liftOver**: https://hgdownload.soe.ucsc.edu/admin/exe/
- **hg19-to-hg38 chain file**: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
- **GRCh38 reference FASTA**: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
- **dbSNP VCF**: https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/
- **Ensembl REST API**: https://rest.ensembl.org/
- **SBayesRC**: https://github.com/zhilizheng/SBayesRC
