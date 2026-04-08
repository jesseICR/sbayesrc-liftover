#!/usr/bin/env python3
"""
liftover.py  --  Lift SBayesRC snp.info from hg19 (GRCh37) to hg38 (GRCh38).

Pipeline:
  1. UCSC liftOver           position-based coordinate conversion
  2. FASTA ref check         flag suspects where neither allele matches hg38 ref
  3. dbSNP cross-validation  catch liftOver errors invisible to FASTA, fill gaps
  4. Ensembl cross-validation  validate unverified positions, rescue remaining
  5. Allele annotation       strand complement + ref/alt determination

Usage:
    bash main.sh
    ENSEMBL_FULL=1 bash main.sh             # exhaustive Ensembl validation (~3h, cached)

Input:  snp.info   (SBayesRC tab-delimited, hg19 coordinates)
Output: sbayesrc_hg38.csv              (clean: chrom, pos, ref, alt, rsid)
        sbayesrc_liftover_results.csv  (verbose: all columns, all SNPs)
"""
import gzip
import json
import os
import platform
import subprocess
import sys
import time
import urllib.request
from urllib.error import HTTPError, URLError

import pandas as pd
import pysam

pd.set_option("future.no_silent_downcasting", True)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.join(ROOT, "tools")
BIN = os.path.join(TOOLS, "bin")
TMP = os.path.join(ROOT, "tmp")

LIFTOVER = os.path.join(BIN, "liftOver")
CHAIN = os.path.join(TOOLS, "hg19ToHg38.over.chain.gz")
HG38_FA = os.path.join(TOOLS, "hg38.fa")
HG38_FA_GZ = HG38_FA + ".gz"
DBSNP_TSV = os.path.join(TOOLS, "dbsnp_lookup.tsv")
ENSEMBL_CACHE = os.path.join(TMP, "ensembl_cache.json")

SNPINFO_URL = "https://github.com/jesseICR/sbayesrc-liftover/releases/download/v1.0/snp.info"
CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
HG38_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
DBSNP_URL = "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"
ENSEMBL_URL = "https://rest.ensembl.org/variation/homo_sapiens"

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

# GRCh38.p14 RefSeq accessions -> chromosome number (autosomes only)
REFSEQ_TO_CHROM = {
    "NC_000001.11": "1",  "NC_000002.12": "2",  "NC_000003.12": "3",
    "NC_000004.12": "4",  "NC_000005.10": "5",  "NC_000006.12": "6",
    "NC_000007.14": "7",  "NC_000008.11": "8",  "NC_000009.12": "9",
    "NC_000010.11": "10", "NC_000011.10": "11", "NC_000012.12": "12",
    "NC_000013.11": "13", "NC_000014.9": "14",  "NC_000015.10": "15",
    "NC_000016.10": "16", "NC_000017.11": "17", "NC_000018.10": "18",
    "NC_000019.10": "19", "NC_000020.11": "20", "NC_000021.9": "21",
    "NC_000022.11": "22",
}


# ---------------------------------------------------------------------------
# Setup  (idempotent — each download skips if the file exists)
# ---------------------------------------------------------------------------
def _liftover_platform():
    s, m = platform.system(), platform.machine()
    plat = {
        ("Linux", "x86_64"): "linux.x86_64",
        ("Darwin", "arm64"): "macOSX.arm64",
        ("Darwin", "x86_64"): "macOSX.x86_64",
    }.get((s, m))
    if plat is None:
        sys.exit(f"Unsupported platform: {s}-{m}")
    return plat


def setup():
    os.makedirs(BIN, exist_ok=True)
    os.makedirs(TMP, exist_ok=True)

    # liftOver binary
    if os.path.isfile(LIFTOVER) and os.access(LIFTOVER, os.X_OK):
        print("  [skip] liftOver binary")
    else:
        print("  [download] liftOver binary ...", flush=True)
        plat = _liftover_platform()
        urllib.request.urlretrieve(
            f"https://hgdownload.soe.ucsc.edu/admin/exe/{plat}/liftOver", LIFTOVER,
        )
        os.chmod(LIFTOVER, 0o755)

    # Chain file
    if os.path.isfile(CHAIN):
        print("  [skip] hg19-to-hg38 chain file")
    else:
        print("  [download] hg19-to-hg38 chain file ...", flush=True)
        urllib.request.urlretrieve(CHAIN_URL, CHAIN)

    # hg38 FASTA (pysam creates .fai automatically on first open)
    if os.path.isfile(HG38_FA):
        print("  [skip] hg38 FASTA")
    else:
        if not os.path.isfile(HG38_FA_GZ):
            print("  [download] hg38.fa.gz (~900 MB) ...", flush=True)
            urllib.request.urlretrieve(HG38_URL, HG38_FA_GZ)
        print("  [decompress] hg38.fa.gz ...", flush=True)
        subprocess.run(["gunzip", "-k", HG38_FA_GZ], check=True)

def build_dbsnp_lookup(rsid_set):
    """Stream dbSNP VCF directly from NCBI FTP, extracting rows that match our rsIDs.
    The full 28 GB VCF is streamed and decompressed on the fly — only the small
    lookup TSV (~200 MB) is saved to disk."""
    if os.path.isfile(DBSNP_TSV):
        print("  [skip] dbSNP lookup table", flush=True)
        return
    print(f"  [build] dbSNP lookup for {len(rsid_set):,} rsIDs "
          f"(streaming ~28 GB VCF, 15-30 min) ...", flush=True)
    # Stream VCF: curl decompresses on the fly via --compressed, but the file
    # is gzipped (not HTTP-compressed), so we pipe through gunzip ourselves.
    proc = subprocess.Popen(
        ["curl", "-sfSL", "--retry", "3", DBSNP_URL],
        stdout=subprocess.PIPE,
    )
    n = 0
    tsv_tmp = DBSNP_TSV + ".tmp"
    with gzip.open(proc.stdout, "rt") as fin, open(tsv_tmp, "w") as fout:
        fout.write("rsid\tchrom\tpos\tref\talt\n")
        for line in fin:
            if line[0] == "#":
                continue
            fields = line.split("\t", 5)
            chrom = REFSEQ_TO_CHROM.get(fields[0])
            if chrom and fields[2] in rsid_set:
                fout.write(f"{fields[2]}\t{chrom}\t{fields[1]}\t{fields[3]}\t{fields[4].split(',')[0]}\n")
                n += 1
    proc.wait()
    if proc.returncode != 0:
        os.remove(tsv_tmp)
        sys.exit("dbSNP VCF download/stream failed")
    os.rename(tsv_tmp, DBSNP_TSV)
    print(f"  [done] {n:,} / {len(rsid_set):,} rsIDs found in dbSNP", flush=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
_FA = None  # lazy-opened pysam.FastaFile


def _open_fasta():
    global _FA
    if _FA is None:
        _FA = pysam.FastaFile(HG38_FA)
    return _FA


def _fasta_ref(df_sub, label=""):
    """Look up the hg38 reference base for rows with chrom + pos_hg38 columns.
    Returns Series of uppercase ref bases indexed by the DataFrame index."""
    fa = _open_fasta()
    chroms = "chr" + df_sub["chrom"].astype(str)
    positions = df_sub["pos_hg38"]
    bases = [
        fa.fetch(c, p - 1, p).upper() for c, p in zip(chroms, positions)
    ]
    return pd.Series(bases, index=df_sub.index)


def _ensembl_post(rsids):
    """POST up to 200 rsIDs to the Ensembl REST API."""
    data = json.dumps({"ids": rsids}).encode()
    req = urllib.request.Request(
        ENSEMBL_URL, data=data,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=30) as resp:
        return json.loads(resp.read().decode())


def _ensembl_query(rsid_list, cache, label=""):
    """Query Ensembl for rsIDs not yet in cache. Saves cache periodically."""
    uncached = [rs for rs in rsid_list if rs not in cache]
    if not uncached:
        return
    n_batches = (len(uncached) + 199) // 200
    print(f"  [query] {len(uncached):,} {label} rsIDs ({n_batches:,} batches) ...",
          flush=True)
    for i in range(0, len(uncached), 200):
        batch = uncached[i : i + 200]
        batch_num = i // 200 + 1
        try:
            result = _ensembl_post(batch)
            for rsid, info in result.items():
                if "error" in info:
                    cache[rsid] = None
                    continue
                hit = None
                for m in info.get("mappings", []):
                    if m.get("assembly_name") == "GRCh38" and m.get("seq_region_name", "").isdigit():
                        hit = {"chrom": int(m["seq_region_name"]), "pos": m["start"]}
                        break
                cache[rsid] = hit
        except (HTTPError, URLError, TimeoutError) as e:
            print(f"    batch {batch_num}: {e}", flush=True)
            for rsid in batch:
                cache.setdefault(rsid, None)
        if batch_num % 500 == 0:
            with open(ENSEMBL_CACHE, "w") as f:
                json.dump(cache, f)
            print(f"    {i + len(batch):,} / {len(uncached):,} "
                  f"({(i + len(batch)) / len(uncached) * 100:.1f}%)", flush=True)
        if i + 200 < len(uncached):
            time.sleep(0.5)
    with open(ENSEMBL_CACHE, "w") as f:
        json.dump(cache, f)


# ---------------------------------------------------------------------------
# Step 1: UCSC liftOver
# ---------------------------------------------------------------------------
def step_liftover(df):
    n = len(df)
    print(f"\n[Step 1] UCSC liftOver ({n:,} SNPs)", flush=True)

    # Write 6-column BED (0-based half-open, with strand for flip detection)
    bed = pd.DataFrame({
        "chrom": "chr" + df["chrom"].astype(str),
        "start": df["pos_hg19"] - 1,
        "end":   df["pos_hg19"],
        "name":  df.index,
        "score": 0,
        "strand": "+",
    })
    bed_in  = os.path.join(TMP, "liftover_in.bed")
    bed_out = os.path.join(TMP, "liftover_out.bed")
    bed_un  = os.path.join(TMP, "liftover_unmapped.bed")
    bed.to_csv(bed_in, sep="\t", header=False, index=False)

    r = subprocess.run(
        [LIFTOVER, bed_in, CHAIN, bed_out, bed_un], capture_output=True, text=True,
    )
    if r.returncode != 0:
        sys.exit(f"liftOver failed: {r.stderr}")

    mapped = pd.read_csv(
        bed_out, sep="\t", header=None,
        names=["chrom38", "start", "end", "idx", "score", "strand_out"],
    )
    mapped["pos_hg38"] = mapped["start"] + 1
    mapped["lo_chrom"] = mapped["chrom38"].str.replace("chr", "", regex=False)
    mapped["strand_flip"] = mapped["strand_out"] == "-"

    # Join back on index (keep first hit if liftOver multi-maps a region)
    mapped = mapped.drop_duplicates("idx").set_index("idx")
    df = df.join(mapped[["pos_hg38", "lo_chrom", "strand_flip"]])

    # Discard chromosome-changed mappings
    has_map = df["pos_hg38"].notna()
    chrom_changed = has_map & (df["chrom"].astype(str) != df["lo_chrom"])
    truly_unmapped = ~has_map

    df.loc[chrom_changed, ["pos_hg38", "strand_flip"]] = [-1, False]
    df["pos_hg38"]    = df["pos_hg38"].fillna(-1).astype(int)
    df["strand_flip"]  = df["strand_flip"].fillna(False).astype(bool)
    df["method"] = ""
    df.loc[df["pos_hg38"] != -1, "method"] = "liftover"
    df.drop(columns=["lo_chrom"], inplace=True)

    n_ok = (df["method"] == "liftover").sum()
    print(f"  Mapped:         {n_ok:>10,}")
    print(f"  Unmapped:       {truly_unmapped.sum():>10,}")
    print(f"  Chrom changed:  {chrom_changed.sum():>10,}")
    print(f"  Strand-flipped: {df['strand_flip'].sum():>10,}")
    return df


# ---------------------------------------------------------------------------
# Step 2: Preliminary FASTA check  (flag liftOver suspects)
# ---------------------------------------------------------------------------
def step_fasta_suspects(df):
    lo = df[df["method"] == "liftover"]
    print(f"\n[Step 2] FASTA ref check ({len(lo):,} liftOver positions)", flush=True)

    ref = _fasta_ref(lo, "suspects")
    a1_match = lo["A1"] == lo.index.map(ref)
    a2_match = lo["A2"] == lo.index.map(ref)
    suspects = set(lo.loc[~a1_match & ~a2_match, "ID"])
    print(f"  Suspects (neither allele matches ref): {len(suspects)}")
    return suspects, ref


# ---------------------------------------------------------------------------
# Step 3: dbSNP cross-validation
# ---------------------------------------------------------------------------
def step_dbsnp(df):
    print(f"\n[Step 3] dbSNP cross-validation", flush=True)

    dbsnp = pd.read_csv(DBSNP_TSV, sep="\t")
    print(f"  Loaded {len(dbsnp):,} dbSNP entries")

    df = df.merge(
        dbsnp[["rsid", "pos"]].drop_duplicates("rsid")
            .rename(columns={"rsid": "ID", "pos": "dbsnp_pos"}),
        on="ID", how="left",
    )
    df["dbsnp_pos"] = df["dbsnp_pos"].fillna(-1).astype(int)

    has_lo = df["pos_hg38"] != -1
    has_db = df["dbsnp_pos"] != -1
    agree    = has_lo & has_db & (df["pos_hg38"] == df["dbsnp_pos"])
    disagree = has_lo & has_db & (df["pos_hg38"] != df["dbsnp_pos"])
    fallback = ~has_lo & has_db

    checked = set(df.loc[has_db, "ID"])

    print(f"  Agree:           {agree.sum():>10,}")
    print(f"  Disagree:        {disagree.sum():>10,}")
    print(f"  Fallback:        {fallback.sum():>10,}")
    print(f"  No dbSNP entry:  {(~has_db).sum():>10,}")

    # Override disagreements — but only if the dbSNP position has a ref allele
    # that matches one of our alleles (or their complements). This guards against
    # cases where dbSNP maps an rsID to a position with a different variant.
    if disagree.any():
        dis = df.loc[disagree].copy()
        # Check FASTA ref at each dbSNP position
        dis_for_ref = pd.DataFrame({
            "chrom": dis["chrom"], "pos_hg38": dis["dbsnp_pos"],
        }, index=dis.index)
        dis_ref = _fasta_ref(dis_for_ref, "dbsnp-check")
        ref_ok = (
            (dis_ref == dis["A1"]) | (dis_ref == dis["A2"]) |
            (dis_ref == dis["A1"].map(COMPLEMENT)) |
            (dis_ref == dis["A2"].map(COMPLEMENT))
        )
        valid = disagree & df.index.isin(dis.index[ref_ok])
        for rsid in df.loc[valid, "ID"]:
            print(f"    override: {rsid}")
        for rsid in df.loc[disagree & ~df.index.isin(dis.index[ref_ok]), "ID"]:
            print(f"    skip (ref mismatch at dbSNP pos): {rsid}")
        df.loc[valid, "pos_hg38"] = df.loc[valid, "dbsnp_pos"]
        df.loc[valid, "method"] = "dbsnp"
        df.loc[valid, "strand_flip"] = False

    # Fill gaps
    df.loc[fallback, "pos_hg38"] = df.loc[fallback, "dbsnp_pos"]
    df.loc[fallback, "method"] = "dbsnp"

    df.drop(columns=["dbsnp_pos"], inplace=True)
    return df, checked


# ---------------------------------------------------------------------------
# Step 4: Ensembl cross-validation
# ---------------------------------------------------------------------------
def step_ensembl(df, suspects, dbsnp_checked):
    still_unmapped      = set(df.loc[df["pos_hg38"] == -1, "ID"])
    unresolved_suspects = suspects & set(df.loc[df["method"] == "liftover", "ID"])
    unverified          = set(df.loc[df["method"] == "liftover", "ID"]) - dbsnp_checked
    high_priority       = still_unmapped | unresolved_suspects

    print(f"\n[Step 4] Ensembl cross-validation", flush=True)
    print(f"  Unmapped:          {len(still_unmapped):>10,}")
    print(f"  Unresolved suspects: {len(unresolved_suspects):>10,}")
    print(f"  Unverified:        {len(unverified):>10,}")

    # Load cache
    cache = {}
    if os.path.isfile(ENSEMBL_CACHE):
        with open(ENSEMBL_CACHE) as f:
            cache = json.load(f)
    print(f"  Cache:             {len(cache):>10,} entries")

    # Always query high-priority
    if high_priority:
        _ensembl_query(list(high_priority), cache, "high-priority")

    # Full validation (opt-in)
    if os.environ.get("ENSEMBL_FULL") == "1":
        _ensembl_query(list(unverified), cache, "unverified")
    else:
        n_uncached = sum(1 for rs in unverified if rs not in cache)
        if n_uncached:
            print(f"  {n_uncached:,} unverified SNPs not in cache "
                  f"(set ENSEMBL_FULL=1 to query, ~3-4 h)")

    # Build lookup from cache -> merge (vectorized, not row-by-row)
    # Iterate cache (small) checking membership, not applicable (huge) checking cache
    applicable = high_priority | unverified
    rows = [(rsid, info["chrom"], info["pos"])
            for rsid, info in cache.items()
            if info is not None and rsid in applicable]

    if not rows:
        print(f"  Nothing to apply")
        n_still = (df["pos_hg38"] == -1).sum()
        print(f"  Unmapped:   {n_still:>6,}", flush=True)
        return df

    ens = pd.DataFrame(rows, columns=["ID", "ens_chrom", "ens_pos"])
    df = df.merge(ens, on="ID", how="left")
    df["ens_pos"]   = df["ens_pos"].fillna(-1).astype(int)
    df["ens_chrom"] = df["ens_chrom"].fillna(-1).astype(int)

    valid    = (df["ens_pos"] != -1) & (df["ens_chrom"] == df["chrom"])
    rescue   = valid & (df["pos_hg38"] == -1)
    override = valid & (df["pos_hg38"] != -1) & (df["pos_hg38"] != df["ens_pos"])

    # Log overrides before applying
    if override.any():
        for _, r in df.loc[override, ["ID", "pos_hg38", "ens_pos"]].iterrows():
            print(f"    override: {r.ID}  "
                  f"liftover={int(r.pos_hg38)}  ensembl={int(r.ens_pos)}")

    df.loc[rescue,   "pos_hg38"] = df.loc[rescue,   "ens_pos"]
    df.loc[rescue,   "method"]   = "ensembl"
    df.loc[override, "pos_hg38"] = df.loc[override, "ens_pos"]
    df.loc[override, "method"]   = "ensembl"
    df.loc[override, "strand_flip"] = False

    df.drop(columns=["ens_pos", "ens_chrom"], inplace=True)

    print(f"  Rescued:    {rescue.sum():>6,}")
    print(f"  Overridden: {override.sum():>6,}")
    print(f"  Unmapped:   {(df['pos_hg38'] == -1).sum():>6,}", flush=True)
    return df


# ---------------------------------------------------------------------------
# Step 5: Final allele annotation
# ---------------------------------------------------------------------------
def step_annotate(df, cached_ref=None):
    print(f"\n[Step 5] Allele annotation", flush=True)

    df["A1_hg38"] = df["A1"]
    df["A2_hg38"] = df["A2"]

    # Complement strand-flipped liftOver alleles.
    # For strand-ambiguous SNPs (A/T, C/G) this swaps labels but not nucleotides —
    # the flip is still real (detected from the chain alignment, not from alleles).
    lo_flip = df["strand_flip"] & (df["method"] == "liftover")
    if lo_flip.any():
        df.loc[lo_flip, "A1_hg38"] = df.loc[lo_flip, "A1"].map(COMPLEMENT)
        df.loc[lo_flip, "A2_hg38"] = df.loc[lo_flip, "A2"].map(COMPLEMENT)

    # Query hg38 FASTA for all mapped positions (reuse step 2 results where possible)
    mapped = df[df["pos_hg38"] != -1]
    if cached_ref is not None:
        # Exclude overridden SNPs — their positions changed since step 2
        still_lo = mapped.index[mapped["method"] == "liftover"]
        valid_cache = cached_ref.loc[cached_ref.index.intersection(still_lo)]
        new = mapped.loc[mapped.index.difference(valid_cache.index)]
        print(f"  FASTA query for {len(new):,} new positions "
              f"({len(valid_cache):,} cached from step 2) ...", flush=True)
        ref = pd.concat([valid_cache, _fasta_ref(new, "final")]) if len(new) else valid_cache
    else:
        print(f"  FASTA query for {len(mapped):,} positions ...", flush=True)
        ref = _fasta_ref(mapped, "final")
    df["ref_hg38"] = df.index.map(ref).fillna("")

    # Non-liftOver SNPs: complement if needed (alleles may be on opposite strand)
    non_lo = df["method"].isin(["dbsnp", "ensembl"])
    if non_lo.any():
        sub = df.loc[non_lo]
        no_match = (sub["A1"] != sub["ref_hg38"]) & (sub["A2"] != sub["ref_hg38"])
        comp_ok  = (sub["A1"].map(COMPLEMENT) == sub["ref_hg38"]) | \
                   (sub["A2"].map(COMPLEMENT) == sub["ref_hg38"])
        needs = no_match & comp_ok
        idx = sub.index[needs]
        if len(idx):
            df.loc[idx, "A1_hg38"] = df.loc[idx, "A1"].map(COMPLEMENT)
            df.loc[idx, "A2_hg38"] = df.loc[idx, "A2"].map(COMPLEMENT)
            df.loc[idx, "strand_flip"] = True
            print(f"  Complemented {len(idx)} non-liftOver SNPs")

    # Determine which allele matches ref
    df["ref_match"] = ""
    has_ref = df["ref_hg38"] != ""
    df.loc[has_ref & (df["A2_hg38"] == df["ref_hg38"]), "ref_match"] = "A2_hg38"
    df.loc[has_ref & (df["A1_hg38"] == df["ref_hg38"]), "ref_match"] = "A1_hg38"
    neither = has_ref & (df["A1_hg38"] != df["ref_hg38"]) & (df["A2_hg38"] != df["ref_hg38"])
    df.loc[neither, "ref_match"] = "neither"

    for val, cnt in df.loc[has_ref, "ref_match"].value_counts().items():
        print(f"    {val:10s}  {cnt:>10,}")
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
RESULTS_COLS = [
    "chrom", "ID", "pos_hg19", "pos_hg38",
    "A1", "A2", "A1_hg38", "A2_hg38",
    "ref_hg38", "ref_match", "strand_flip", "method",
    "Index", "GenPos", "A1Freq", "N", "Block",
]


def main():
    input_file = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "snp.info")

    print("=" * 60)
    print("SBayesRC snp.info  hg19 -> hg38 liftover")
    print("=" * 60, flush=True)

    # Setup
    setup()

    # Download snp.info from GitHub Release if not present
    if not os.path.isfile(input_file):
        print(f"\n  [download] {os.path.basename(input_file)} (~381 MB) ...", flush=True)
        urllib.request.urlretrieve(SNPINFO_URL, input_file)

    print(f"\nReading {os.path.basename(input_file)} ...", flush=True)
    df = pd.read_csv(input_file, sep="\t")
    df = df.rename(columns={"Chrom": "chrom", "PhysPos": "pos_hg19"})
    print(f"  {len(df):,} SNPs on {df['chrom'].nunique()} chromosomes")
    build_dbsnp_lookup(set(df["ID"]))

    # Pipeline
    df = step_liftover(df)                                  # 1
    suspects, cached_ref = step_fasta_suspects(df)          # 2
    df, dbsnp_checked = step_dbsnp(df)                      # 3
    df = step_ensembl(df, suspects, dbsnp_checked)          # 4
    df = step_annotate(df, cached_ref)                      # 5

    # Summary
    n_un = (df["pos_hg38"] == -1).sum()
    print(f"\n{'=' * 60}")
    print(f"Result: {len(df):,} SNPs")
    for m, c in df.loc[df["method"] != "", "method"].value_counts().items():
        print(f"  {m:12s} {c:>10,}")
    print(f"  {'unmapped':12s} {n_un:>10,}")
    print(f"{'=' * 60}", flush=True)

    # Write verbose liftover results (all SNPs, all columns)
    results_csv = os.path.join(ROOT, "sbayesrc_liftover_results.csv")
    df[RESULTS_COLS].to_csv(results_csv, index=False)
    print(f"\nWritten to {os.path.basename(results_csv)}")

    # Write clean hg38 output (mapped SNPs only: chrom, pos, ref, alt, rsid)
    mapped = df[df["pos_hg38"] != -1].copy()
    # alt = whichever allele is NOT the reference
    alt = mapped["A2_hg38"].where(mapped["ref_match"] == "A1_hg38", mapped["A1_hg38"])
    clean = pd.DataFrame({
        "chrom": mapped["chrom"],
        "pos": mapped["pos_hg38"],
        "ref": mapped["ref_hg38"],
        "alt": alt,
        "rsid": mapped["ID"],
    })
    clean_csv = os.path.join(ROOT, "sbayesrc_hg38.csv")
    clean.to_csv(clean_csv, index=False)
    print(f"Written to {os.path.basename(clean_csv)} ({len(clean):,} SNPs)")


if __name__ == "__main__":
    main()
