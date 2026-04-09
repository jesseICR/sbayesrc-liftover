#!/usr/bin/env python3
"""
liftover.py  --  Lift SBayesRC snp.info from hg19 (GRCh37) to hg38 (GRCh38).

Pipeline:
  1. UCSC liftOver           position-based coordinate conversion (hg19 -> hg38)
  2. dbSNP cross-validation  confirm liftOver positions, rescue failures, discard conflicts
  3. FASTA reference check   verify ref/alt alleles against hg38 FASTA and dbSNP
  4. Allele annotation       strand complement + ref/alt determination
  5. 1000G EUR validation    compare allele frequencies against 1000 Genomes (informational)

Inclusion criteria -- every SNP in sbayesrc_hg38.csv must satisfy ALL of:
  1. Its rsID exists in dbSNP with a matching hg38 chrom and position
  2. If UCSC liftOver succeeded, the liftOver position must equal the dbSNP position
  3. The dbSNP reference allele must match the hg38 FASTA reference base
  4. At least one allele (A1/A2, accounting for strand) must match the FASTA ref
  5. The other (non-ref) allele must appear in dbSNP's alt allele(s) for that rsID

Usage:
    bash main.sh

Input:  snp.info   (SBayesRC tab-delimited, hg19 coordinates)
Output: sbayesrc_hg38.csv              (clean: chrom, pos, ref, alt, rsid)
        sbayesrc_liftover_results.csv  (verbose: all columns, all SNPs)
        kg_validation/allele_freq_validation.png  (1000G frequency scatter plot)
"""
import gzip
import io
import os
import platform
import subprocess
import sys
import urllib.request

import pandas as pd
import pysam

pd.set_option("future.no_silent_downcasting", True)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = os.getcwd()
TOOLS = os.path.join(ROOT, "tools")
BIN = os.path.join(TOOLS, "bin")
TMP = os.path.join(ROOT, "tmp")

LIFTOVER = os.path.join(BIN, "liftOver")
CHAIN = os.path.join(TOOLS, "hg19ToHg38.over.chain.gz")
HG38_FA = os.path.join(TOOLS, "hg38.fa")
HG38_FA_GZ = HG38_FA + ".gz"
DBSNP_TSV = os.path.join(TOOLS, "dbsnp_lookup.tsv")
KG_PVAR_ZST = os.path.join(TOOLS, "kg_all.pvar.zst")
KG_LOOKUP_TSV = os.path.join(TOOLS, "kg_eur_lookup.tsv")

SNPINFO_URL = "https://github.com/jesseICR/sbayesrc-liftover/releases/download/v1.0/snp.info"
CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
HG38_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
DBSNP_URL = "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"
KG_PVAR_URL = "https://www.dropbox.com/scl/fi/fn0bcm5oseyuawxfvkcpb/all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz&dl=1"

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
# Setup  (idempotent -- each download skips if the file exists)
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

    # 1000 Genomes pvar (contains pre-computed AF_EUR_unrel in INFO)
    if os.path.isfile(KG_PVAR_ZST):
        print("  [skip] 1000G pvar")
    else:
        print("  [download] 1000G pvar.zst (~700 MB) ...", flush=True)
        subprocess.run(
            ["curl", "-fSL", "-o", KG_PVAR_ZST, KG_PVAR_URL], check=True,
        )

def build_dbsnp_lookup(rsid_set):
    """Stream dbSNP VCF directly from NCBI FTP, extracting rows that match our rsIDs.
    The full 28 GB VCF is streamed and decompressed on the fly -- only the small
    lookup TSV (~200 MB) is saved to disk.

    Stores the full comma-separated ALT field (e.g. "A,G" for multi-allelic sites)
    so that the alt allele check in step 3 can verify against all possible alts."""
    if os.path.isfile(DBSNP_TSV):
        print("  [skip] dbSNP lookup table", flush=True)
        return
    print(f"  [build] dbSNP lookup for {len(rsid_set):,} rsIDs "
          f"(streaming ~28 GB VCF, 15-30 min) ...", flush=True)
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
                fout.write(f"{fields[2]}\t{chrom}\t{fields[1]}\t{fields[3]}\t{fields[4]}\n")
                n += 1
    proc.wait()
    if proc.returncode != 0:
        os.remove(tsv_tmp)
        sys.exit("dbSNP VCF download/stream failed")
    os.rename(tsv_tmp, DBSNP_TSV)
    print(f"  [done] {n:,} / {len(rsid_set):,} rsIDs found in dbSNP", flush=True)


def build_kg_lookup(rsid_set):
    """Stream 1000G pvar.zst and extract AF_EUR_unrel for matching rsIDs."""
    if os.path.isfile(KG_LOOKUP_TSV):
        print("  [skip] 1000G EUR lookup table", flush=True)
        return
    print(f"  [build] 1000G EUR lookup for {len(rsid_set):,} rsIDs "
          f"(streaming pvar.zst) ...", flush=True)
    import zstandard
    dctx = zstandard.ZstdDecompressor()
    n = 0
    tsv_tmp = KG_LOOKUP_TSV + ".tmp"
    with open(KG_PVAR_ZST, "rb") as fin, open(tsv_tmp, "w") as fout:
        fout.write("rsid\tchrom\tpos\tref\talt\taf_eur_unrel\n")
        reader = dctx.stream_reader(fin)
        text = io.TextIOWrapper(reader, encoding="utf-8")
        for line in text:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            rsid = fields[2]
            if rsid not in rsid_set:
                continue
            info = fields[7]
            af = ""
            for part in info.split(";"):
                if part.startswith("AF_EUR_unrel="):
                    af = part.split("=", 1)[1]
                    break
            if af:
                fout.write(f"{rsid}\t{fields[0]}\t{fields[1]}\t{fields[3]}\t{fields[4]}\t{af}\n")
                n += 1
    os.rename(tsv_tmp, KG_LOOKUP_TSV)
    print(f"  [done] {n:,} / {len(rsid_set):,} rsIDs found in 1000G EUR", flush=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
_FA = None  # lazy-opened pysam.FastaFile


def _open_fasta():
    global _FA
    if _FA is None:
        _FA = pysam.FastaFile(HG38_FA)
    return _FA


def _fasta_ref(df_sub):
    """Look up the hg38 reference base for rows with chrom + pos_hg38 columns.
    Returns Series of uppercase ref bases indexed by the DataFrame index."""
    fa = _open_fasta()
    chroms = "chr" + df_sub["chrom"].astype(str)
    positions = df_sub["pos_hg38"]
    bases = [
        fa.fetch(c, p - 1, p).upper() for c, p in zip(chroms, positions)
    ]
    return pd.Series(bases, index=df_sub.index)


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
    mapped["pos_liftover"] = mapped["start"] + 1
    mapped["lo_chrom"] = mapped["chrom38"].str.replace("chr", "", regex=False)
    mapped["strand_flip"] = mapped["strand_out"] == "-"

    # Join back on index (keep first hit if liftOver multi-maps a region)
    mapped = mapped.drop_duplicates("idx").set_index("idx")
    df = df.join(mapped[["pos_liftover", "lo_chrom", "strand_flip"]])

    # Discard chromosome-changed mappings
    has_map = df["pos_liftover"].notna()
    chrom_changed = has_map & (df["chrom"].astype(str) != df["lo_chrom"])
    truly_unmapped = ~has_map

    df.loc[chrom_changed, ["pos_liftover", "strand_flip"]] = [-1, False]
    df["pos_liftover"] = df["pos_liftover"].fillna(-1).astype(int)
    df["strand_flip"]  = df["strand_flip"].fillna(False).astype(bool)
    df.drop(columns=["lo_chrom"], inplace=True)

    n_ok = (df["pos_liftover"] != -1).sum()
    n_flip = df["strand_flip"].sum()
    print(f"  Mapped:         {n_ok:>10,}")
    print(f"  Unmapped:       {truly_unmapped.sum():>10,}")
    print(f"  Chrom changed:  {chrom_changed.sum():>10,}")
    print(f"  Strand-flipped: {n_flip:>10,}")
    return df


# ---------------------------------------------------------------------------
# Step 2: dbSNP cross-validation
# ---------------------------------------------------------------------------
def step_dbsnp(df):
    print(f"\n[Step 2] dbSNP cross-validation", flush=True)

    dbsnp = pd.read_csv(DBSNP_TSV, sep="\t", dtype={"chrom": str})
    print(f"  Loaded {len(dbsnp):,} dbSNP entries")

    # Merge dbSNP chrom, pos, ref, alt onto the dataframe
    dbsnp_dedup = (
        dbsnp[["rsid", "chrom", "pos", "ref", "alt"]]
        .drop_duplicates("rsid")
        .rename(columns={
            "rsid": "ID",
            "chrom": "dbsnp_chrom",
            "pos": "pos_dbsnp",
            "ref": "dbsnp_ref",
            "alt": "dbsnp_alt",
        })
    )
    df = df.merge(dbsnp_dedup, on="ID", how="left")
    df["pos_dbsnp"]    = df["pos_dbsnp"].fillna(-1).astype(int)
    df["dbsnp_chrom"]  = df["dbsnp_chrom"].fillna("")
    df["dbsnp_ref"]    = df["dbsnp_ref"].fillna("")
    df["dbsnp_alt"]    = df["dbsnp_alt"].fillna("")

    # Classify each SNP (mutually exclusive, exhaustive)
    has_lo = df["pos_liftover"] != -1
    has_db = df["pos_dbsnp"] != -1
    db_chrom_ok = df["chrom"].astype(str) == df["dbsnp_chrom"]
    db_usable = has_db & db_chrom_ok

    confirmed = has_lo & db_usable & (df["pos_liftover"] == df["pos_dbsnp"])
    conflict  = has_lo & has_db & ~confirmed   # liftOver mapped but dbSNP disagrees (pos or chrom)
    rescue    = ~has_lo & db_usable
    no_dbsnp  = has_lo & ~has_db               # liftOver worked but no dbSNP to confirm
    unmapped  = ~has_lo & ~db_usable           # neither source provides a usable position

    df["status"] = ""
    df.loc[confirmed, "status"] = "confirmed"
    df.loc[conflict,  "status"] = "conflict"
    df.loc[rescue,    "status"] = "rescue"
    df.loc[no_dbsnp,  "status"] = "no_dbsnp"
    df.loc[unmapped,  "status"] = "unmapped"

    # Set pos_hg38: confirmed use liftOver pos, rescue use dbSNP pos, rest get -1
    df["pos_hg38"] = -1
    df.loc[confirmed, "pos_hg38"] = df.loc[confirmed, "pos_liftover"]
    df.loc[rescue,    "pos_hg38"] = df.loc[rescue,    "pos_dbsnp"]

    print(f"  Confirmed (liftOver = dbSNP):              {confirmed.sum():>10,}")
    print(f"  Conflict  (liftOver != dbSNP):             {conflict.sum():>10,}")
    print(f"  Rescue    (liftOver failed, dbSNP rescue): {rescue.sum():>10,}")
    print(f"  No dbSNP  (liftOver OK, no dbSNP entry):  {no_dbsnp.sum():>10,}")
    print(f"  Unmapped  (neither source):                {unmapped.sum():>10,}")

    df.drop(columns=["dbsnp_chrom"], inplace=True)
    return df


# ---------------------------------------------------------------------------
# Step 3: FASTA reference validation
# ---------------------------------------------------------------------------
def step_fasta_validation(df):
    checkable = df["status"].isin(["confirmed", "rescue"])
    n_check = checkable.sum()
    print(f"\n[Step 3] FASTA reference validation ({n_check:,} SNPs)", flush=True)

    # Fetch the hg38 FASTA ref base at each position
    sub = df.loc[checkable]
    fasta_ref = _fasta_ref(sub)
    df["fasta_ref"] = ""
    df.loc[checkable, "fasta_ref"] = fasta_ref

    # Check 1: dbSNP ref must match FASTA ref
    fasta_mismatch = checkable & (df["dbsnp_ref"] != df["fasta_ref"])
    n_fasta_mismatch = fasta_mismatch.sum()
    df.loc[fasta_mismatch, "status"] = "fasta_mismatch"
    df.loc[fasta_mismatch, "pos_hg38"] = -1

    # Check 2: at least one allele (or complement) matches FASTA ref
    still_ok = df["status"].isin(["confirmed", "rescue"])
    allele_ok = (
        (df["A1"] == df["fasta_ref"]) | (df["A2"] == df["fasta_ref"]) |
        (df["A1"].map(COMPLEMENT) == df["fasta_ref"]) |
        (df["A2"].map(COMPLEMENT) == df["fasta_ref"])
    )
    allele_mismatch = still_ok & ~allele_ok
    n_allele_mismatch = allele_mismatch.sum()
    df.loc[allele_mismatch, "status"] = "allele_mismatch"
    df.loc[allele_mismatch, "pos_hg38"] = -1

    # Check 3: the non-ref allele must appear in dbSNP's alt allele(s).
    #
    # First we figure out which strand we're on by seeing which allele
    # (A1/A2 or their complement) matched the FASTA ref in check 2.
    # The OTHER allele on that same strand is the non-ref allele.
    # That non-ref allele must exist in dbSNP's comma-separated alt list.
    #
    # Example (forward strand): A1=C matches ref, so A2=T is the non-ref.
    #   Check T is in dbSNP alts.
    # Example (reverse strand): complement(A1)=G matches ref, so
    #   complement(A2)=A is the non-ref. Check A is in dbSNP alts.
    still_ok = df["status"].isin(["confirmed", "rescue"])
    dbsnp_alt_sets = df["dbsnp_alt"].str.split(",")

    def _alt_in_dbsnp(row):
        if not row["_still_ok"]:
            return True  # skip already-excluded rows
        alts = row["_dbsnp_alts"]
        if not isinstance(alts, list):
            return False
        ref = row["fasta_ref"]
        a1, a2 = row["A1"], row["A2"]
        a1c = COMPLEMENT.get(a1, "")
        a2c = COMPLEMENT.get(a2, "")
        # Identify which strand, then pick the non-ref allele on that strand
        if a1 == ref:
            non_ref = a2          # forward strand: A2 is alt
        elif a2 == ref:
            non_ref = a1          # forward strand: A1 is alt
        elif a1c == ref:
            non_ref = a2c         # reverse strand: complement(A2) is alt
        elif a2c == ref:
            non_ref = a1c         # reverse strand: complement(A1) is alt
        else:
            return False          # shouldn't happen -- check 2 already verified
        return non_ref in alts

    check_df = pd.DataFrame({
        "A1": df["A1"], "A2": df["A2"],
        "fasta_ref": df["fasta_ref"],
        "_dbsnp_alts": dbsnp_alt_sets,
        "_still_ok": still_ok,
    })
    alt_ok = check_df.apply(_alt_in_dbsnp, axis=1)
    alt_mismatch = still_ok & ~alt_ok
    n_alt_mismatch = alt_mismatch.sum()
    df.loc[alt_mismatch, "status"] = "alt_mismatch"
    df.loc[alt_mismatch, "pos_hg38"] = -1

    n_passed = df["status"].isin(["confirmed", "rescue"]).sum()
    print(f"  dbSNP ref matches FASTA ref:  {(n_check - n_fasta_mismatch):>10,}")
    print(f"  dbSNP ref != FASTA ref:       {n_fasta_mismatch:>10,}")
    print(f"  Allele matches FASTA ref:     {(n_check - n_fasta_mismatch - n_allele_mismatch):>10,}")
    print(f"  Neither allele matches ref:   {n_allele_mismatch:>10,}")
    print(f"  Alt allele in dbSNP alts:     {(n_check - n_fasta_mismatch - n_allele_mismatch - n_alt_mismatch):>10,}")
    print(f"  Alt allele not in dbSNP:      {n_alt_mismatch:>10,}")
    print(f"  Passed all validation:        {n_passed:>10,}")
    return df


# ---------------------------------------------------------------------------
# Step 4: Allele annotation
# ---------------------------------------------------------------------------
def step_annotate(df):
    print(f"\n[Step 4] Allele annotation", flush=True)

    passed = df["status"].isin(["confirmed", "rescue"])
    df["A1_hg38"] = df["A1"]
    df["A2_hg38"] = df["A2"]

    # Complement strand-flipped liftOver alleles (confirmed SNPs only)
    lo_flip = passed & (df["status"] == "confirmed") & df["strand_flip"]
    if lo_flip.any():
        df.loc[lo_flip, "A1_hg38"] = df.loc[lo_flip, "A1"].map(COMPLEMENT)
        df.loc[lo_flip, "A2_hg38"] = df.loc[lo_flip, "A2"].map(COMPLEMENT)
        print(f"  Complemented {lo_flip.sum()} strand-flipped liftOver SNPs")

    # Rescue SNPs: complement if needed (alleles may be on opposite strand)
    rescue = passed & (df["status"] == "rescue")
    if rescue.any():
        sub = df.loc[rescue]
        ref = df.loc[rescue, "fasta_ref"]
        no_match = (sub["A1"] != ref) & (sub["A2"] != ref)
        comp_ok  = (sub["A1"].map(COMPLEMENT) == ref) | \
                   (sub["A2"].map(COMPLEMENT) == ref)
        needs = rescue & df.index.isin(sub.index[no_match & comp_ok])
        if needs.any():
            df.loc[needs, "A1_hg38"] = df.loc[needs, "A1"].map(COMPLEMENT)
            df.loc[needs, "A2_hg38"] = df.loc[needs, "A2"].map(COMPLEMENT)
            df.loc[needs, "strand_flip"] = True
            print(f"  Complemented {needs.sum()} rescue SNPs")

    # Determine which allele matches ref
    df["ref_match"] = ""
    df.loc[passed & (df["A1_hg38"] == df["fasta_ref"]), "ref_match"] = "A1_hg38"
    df.loc[passed & (df["A2_hg38"] == df["fasta_ref"]), "ref_match"] = "A2_hg38"

    for val, cnt in df.loc[passed, "ref_match"].value_counts().items():
        print(f"    {val:10s}  {cnt:>10,}")
    return df


# ---------------------------------------------------------------------------
# Step 5: 1000G EUR allele frequency validation
# ---------------------------------------------------------------------------
def step_kg_validation(df):
    """Compare allele frequencies against 1000 Genomes EUR as a sanity check.

    Both our data (A1_hg38/A2_hg38) and the 1000G pvar are on the hg38
    forward strand, so we match by exact allele identity -- no strand
    complement logic is needed. For multi-allelic sites in 1000G (multiple
    rows per rsID), this naturally picks the row with the matching alt.
    """
    print(f"\n[Step 5] 1000G EUR allele frequency validation", flush=True)

    kg = pd.read_csv(KG_LOOKUP_TSV, sep="\t", dtype={"chrom": str, "pos": int})
    print(f"  Loaded {len(kg):,} 1000G entries")

    passed = df["status"].isin(["confirmed", "rescue"])

    # Determine our alt allele (the non-reference allele on hg38+ strand).
    # ref_match == "A1_hg38" means A1 is ref, so A2 is alt; and vice versa.
    df["_snp_alt"] = ""
    is_a1_ref = passed & (df["ref_match"] == "A1_hg38")
    is_a2_ref = passed & (df["ref_match"] == "A2_hg38")
    df.loc[is_a1_ref, "_snp_alt"] = df.loc[is_a1_ref, "A2_hg38"]
    df.loc[is_a2_ref, "_snp_alt"] = df.loc[is_a2_ref, "A1_hg38"]

    # Merge on all five columns: rsID, chrom, pos, ref, alt.
    # For multi-allelic 1000G sites this picks the row whose alt matches
    # our alt allele, rather than arbitrarily taking the first row.
    # kg_af is the frequency of kg_alt in 1000G EUR unrelated.
    kg_merge = kg.rename(columns={
        "rsid": "ID",
        "pos": "pos_hg38",
        "ref": "fasta_ref",
        "alt": "_snp_alt",
        "af_eur_unrel": "kg_af",
    })[["ID", "chrom", "pos_hg38", "fasta_ref", "_snp_alt", "kg_af"]]
    kg_merge["chrom"] = kg_merge["chrom"].astype(int)
    kg_merge["pos_hg38"] = kg_merge["pos_hg38"].astype(int)
    df = df.merge(kg_merge, on=["ID", "chrom", "pos_hg38", "fasta_ref", "_snp_alt"],
                  how="left")

    in_kg = passed & df["kg_af"].notna()
    not_in_kg = passed & df["kg_af"].isna()

    # Break down why SNPs weren't matched in 1000G
    kg_rsids = set(kg["rsid"])
    not_in_kg_no_rsid = not_in_kg & ~df["ID"].isin(kg_rsids)
    not_in_kg_no_allele = not_in_kg & df["ID"].isin(kg_rsids)

    print(f"  Matched in 1000G:             {in_kg.sum():>10,}")
    print(f"  Not in 1000G (rsID missing):  {not_in_kg_no_rsid.sum():>10,}")
    print(f"  Not in 1000G (allele mismatch): {not_in_kg_no_allele.sum():>10,}")

    # Compute a1_freq_kg: frequency of A1 (effect allele) in 1000G EUR.
    # kg_af is the frequency of the alt allele. So:
    #   A1 is our alt (ref_match == "A2_hg38") → a1_freq_kg = kg_af
    #   A1 is our ref (ref_match == "A1_hg38") → a1_freq_kg = 1 - kg_af
    df["a1_freq_kg"] = float("nan")
    a1_is_alt = in_kg & (df["ref_match"] == "A2_hg38")
    a1_is_ref = in_kg & (df["ref_match"] == "A1_hg38")
    df.loc[a1_is_alt, "a1_freq_kg"] = df.loc[a1_is_alt, "kg_af"]
    df.loc[a1_is_ref, "a1_freq_kg"] = 1 - df.loc[a1_is_ref, "kg_af"]

    print(f"  Allele-frequency assigned:    {in_kg.sum():>10,}")

    # Flag large frequency differences (|A1Freq - a1_freq_kg| > 0.2)
    has_freq = df["a1_freq_kg"].notna()
    freq_diff = (df.loc[has_freq, "A1Freq"] - df.loc[has_freq, "a1_freq_kg"]).abs()
    large_diff = freq_diff > 0.2
    n_large = large_diff.sum()
    print(f"  |freq diff| > 0.2:            {n_large:>10,}")

    if n_large > 0:
        bad_idx = freq_diff.index[large_diff]
        print(f"\n  rsIDs with |A1Freq - a1_freq_kg| > 0.2:")
        for _, r in df.loc[bad_idx, ["ID", "chrom", "pos_hg38", "A1Freq", "a1_freq_kg"]].head(50).iterrows():
            print(f"    {r.ID}  chr{r.chrom}:{int(r.pos_hg38)}  "
                  f"A1Freq={r.A1Freq:.4f}  kg={r.a1_freq_kg:.4f}  "
                  f"diff={abs(r.A1Freq - r.a1_freq_kg):.4f}")
        if n_large > 50:
            print(f"    ... and {n_large - 50} more")

    # Not-in-1000G list (log first 20)
    if not_in_kg.any():
        print(f"\n  Sample rsIDs not matched in 1000G (first 20):")
        for _, r in df.loc[not_in_kg, ["ID", "chrom", "pos_hg38"]].head(20).iterrows():
            print(f"    {r.ID}  chr{r.chrom}:{int(r.pos_hg38)}")

    # Scatter plot (only SNPs with a valid a1_freq_kg; unmatched SNPs
    # are excluded entirely, not plotted at frequency 0)
    kg_dir = os.path.join(ROOT, "kg_validation")
    os.makedirs(kg_dir, exist_ok=True)
    plot_path = os.path.join(kg_dir, "allele_freq_validation.png")
    _make_freq_plot(df.loc[has_freq], plot_path)
    print(f"\n  Scatter plot: kg_validation/{os.path.basename(plot_path)}")

    df.drop(columns=["_snp_alt", "kg_af"], inplace=True)
    return df


def _make_freq_plot(df, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = df["A1Freq"].values
    y = df["a1_freq_kg"].values

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(x, y, s=0.1, alpha=0.05, rasterized=True)
    ax.plot([0, 1], [0, 1], "r--", linewidth=1, alpha=0.5)
    ax.set_xlabel("A1 freq (SBayesRC snp.info)")
    ax.set_ylabel("A1 freq (1000G EUR unrelated)")
    ax.set_title(f"Allele frequency validation (n={len(df):,})")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
RESULTS_COLS = [
    "chrom", "ID", "pos_hg19", "pos_hg38", "pos_liftover", "pos_dbsnp",
    "A1", "A2", "A1_hg38", "A2_hg38",
    "dbsnp_ref", "dbsnp_alt", "fasta_ref", "ref_match", "strand_flip", "status",
    "a1_freq_kg",
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
    rsid_set = set(df["ID"])
    build_dbsnp_lookup(rsid_set)
    build_kg_lookup(rsid_set)

    # Pipeline
    df = step_liftover(df)           # 1
    df = step_dbsnp(df)              # 2
    df = step_fasta_validation(df)   # 3
    df = step_annotate(df)           # 4

    # Check for duplicate chrom+pos among passed SNPs before 1000G validation.
    # Two different rsIDs mapping to the same hg38 position is ambiguous.
    passed = df["status"].isin(["confirmed", "rescue"])
    dup_mask = passed & df.duplicated(subset=["chrom", "pos_hg38"], keep=False)
    n_dup = dup_mask.sum()
    if n_dup > 0:
        df.loc[dup_mask, "status"] = "duplicate_pos"
        df.loc[dup_mask, "pos_hg38"] = -1
        print(f"\n  Duplicate chrom+pos excluded: {n_dup:,}")
        for _, r in df.loc[dup_mask, ["ID", "chrom", "pos_hg38"]].head(20).iterrows():
            print(f"    {r.ID}  chr{r.chrom}:{int(r.pos_hg38)}")
    else:
        print(f"\n  No duplicate chrom+pos among passed SNPs")

    df = step_kg_validation(df)      # 5

    # ---- Final summary -------------------------------------------------------
    n_total = len(df)
    passed = df["status"].isin(["confirmed", "rescue"])
    n_confirmed = (df["status"] == "confirmed").sum()
    n_rescue    = (df["status"] == "rescue").sum()
    n_conflict  = (df["status"] == "conflict").sum()
    n_no_dbsnp  = (df["status"] == "no_dbsnp").sum()
    n_unmapped  = (df["status"] == "unmapped").sum()
    n_fasta     = (df["status"] == "fasta_mismatch").sum()
    n_allele    = (df["status"] == "allele_mismatch").sum()
    n_alt       = (df["status"] == "alt_mismatch").sum()
    n_dup       = (df["status"] == "duplicate_pos").sum()

    print(f"\n{'=' * 60}")
    print(f"FINAL SUMMARY: {n_total:,} SNPs")
    print(f"{'=' * 60}")
    print(f"\n  Included in sbayesrc_hg38.csv:")
    print(f"    confirmed (liftOver + dbSNP agree):  {n_confirmed:>10,}")
    print(f"    rescue    (dbSNP only):              {n_rescue:>10,}")
    print(f"    TOTAL INCLUDED:                      {n_confirmed + n_rescue:>10,}")
    n_excluded = n_conflict + n_no_dbsnp + n_unmapped + n_fasta + n_allele + n_alt + n_dup
    print(f"\n  Excluded:")
    print(f"    conflict  (liftOver != dbSNP):       {n_conflict:>10,}")
    print(f"    no_dbsnp  (rsID not in dbSNP):       {n_no_dbsnp:>10,}")
    print(f"    unmapped  (neither source):          {n_unmapped:>10,}")
    print(f"    fasta_mismatch (dbSNP ref != FASTA): {n_fasta:>10,}")
    print(f"    allele_mismatch (no allele = ref):   {n_allele:>10,}")
    print(f"    alt_mismatch (alt not in dbSNP):     {n_alt:>10,}")
    print(f"    duplicate_pos (same chrom+pos):      {n_dup:>10,}")
    print(f"    TOTAL EXCLUDED:                      {n_excluded:>10,}")
    print(f"{'=' * 60}", flush=True)

    # ---- Write verbose liftover results (all SNPs, all columns) ---------------
    results_csv = os.path.join(ROOT, "sbayesrc_liftover_results.csv")
    df[RESULTS_COLS].to_csv(results_csv, index=False)
    print(f"\nWritten to {os.path.basename(results_csv)}")

    # ---- Write clean hg38 output (passed SNPs only) --------------------------
    mapped = df[passed].copy()
    # alt = whichever allele is NOT the reference
    alt = mapped["A2_hg38"].where(mapped["ref_match"] == "A1_hg38", mapped["A1_hg38"])
    clean = pd.DataFrame({
        "chrom": mapped["chrom"],
        "pos": mapped["pos_hg38"],
        "ref": mapped["fasta_ref"],
        "alt": alt,
        "rsid": mapped["ID"],
    })
    clean_csv = os.path.join(ROOT, "sbayesrc_hg38.csv")
    clean.to_csv(clean_csv, index=False)
    print(f"Written to {os.path.basename(clean_csv)} ({len(clean):,} SNPs)")


if __name__ == "__main__":
    main()
