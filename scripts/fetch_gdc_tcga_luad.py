#!/usr/bin/env python
"""
fetch_gdc_tcga_luad.py — build a uniformly processed TCGA-LUAD tumor vs adjacent-normal dataset
================================================================================================
Downloads the GDC "STAR - Counts" gene-expression files of the TCGA-LUAD project through the
public GDC REST API (no account or token needed: these files are open access) and assembles

    data/raw/tcga_gdc/counts_matrix.csv   genes x samples, raw unstranded STAR counts
    data/raw/tcga_gdc/tpm_matrix.csv      genes x samples, STAR tpm_unstranded
    data/raw/tcga_gdc/sample_sheet.csv    sample_id, patient_id, group (tumor | normal), sample_type, ...
    data/raw/tcga_gdc/gdc_files_manifest.csv   the file list returned by GDC (provenance)

Both groups come from the SAME project processed by the SAME GDC pipeline (STAR alignment, same
gene model), which removes the cohort / quantification confounding of a TCGA-vs-GTEx comparison.

GDC query (POST https://api.gdc.cancer.gov/files)
    cases.project.project_id          = TCGA-LUAD
    data_type                          = Gene Expression Quantification
    analysis.workflow_type             = STAR - Counts
    cases.samples.sample_type          in {Primary Tumor, Solid Tissue Normal}
Each file is fetched with GET https://api.gdc.cancer.gov/data/<file_id>.  A STAR counts file is a
tab-separated table (about 4.2 MB) with one comment line, a header
    gene_id gene_name gene_type unstranded stranded_first stranded_second
    tpm_unstranded fpkm_unstranded fpkm_uq_unstranded
and four summary rows (N_unmapped, N_multimapping, N_noFeature, N_ambiguous) that are skipped.

Usage
    python scripts/fetch_gdc_tcga_luad.py                       # everything (~600 files, ~2.5 GB)
    python scripts/fetch_gdc_tcga_luad.py --paired-only         # normals + their matched tumors (~0.5 GB)
    python scripts/fetch_gdc_tcga_luad.py --workers 8 --out data/raw/tcga_gdc
    python scripts/fetch_gdc_tcga_luad.py --all-genes           # keep non-coding genes too

Then set in config.py:  DATA_SOURCE = "tcga_gdc"  and  DE_METHOD = "paired" (or "limma_voom").
See docs/tcga_paired_design.md.  Downloaded files are cached in <out>/files/ so the script can be
re-run safely; use --force to download again.
"""
import argparse
import json
import logging
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

GDC_FILES_URL = "https://api.gdc.cancer.gov/files"
GDC_DATA_URL = "https://api.gdc.cancer.gov/data/"
SAMPLE_TYPES = {"Primary Tumor": "tumor", "Solid Tissue Normal": "normal"}
USER_AGENT = "gloom-fetch-gdc/0.2"

log = logging.getLogger("fetch_gdc_tcga_luad")


# ── GDC query ──────────────────────────────────────────────────────────────────────────────────────

def _post_json(url, payload, timeout=60, retries=3):
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, data=json.dumps(payload).encode("utf-8"),
                headers={"Content-Type": "application/json", "Accept": "application/json",
                         "User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except Exception as exc:
            last = exc
            time.sleep(1.5 ** attempt)
    raise RuntimeError(f"GDC query failed after {retries} attempts: {last}")


def query_files(project="TCGA-LUAD"):
    """Return the GDC file records for STAR counts of tumor / adjacent-normal samples."""
    filters = {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
        {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
        {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
        {"op": "in", "content": {"field": "cases.samples.sample_type", "value": list(SAMPLE_TYPES)}},
    ]}
    payload = {
        "filters": filters,
        "fields": ",".join([
            "file_id", "file_name", "file_size",
            "cases.submitter_id",
            "cases.samples.submitter_id", "cases.samples.sample_type",
        ]),
        "format": "JSON",
        "size": "2000",
    }
    resp = _post_json(GDC_FILES_URL, payload)
    hits = resp["data"]["hits"]
    total = resp["data"]["pagination"]["total"]
    if total > len(hits):
        log.warning(f"GDC reports {total} files but only {len(hits)} were returned.")
    rows = []
    for h in hits:
        found = None                           # first (case, sample) of a wanted sample type
        for case in h.get("cases", []):
            for samp in case.get("samples", []):
                if samp.get("sample_type") in SAMPLE_TYPES:
                    found = (case, samp)
                    break
            if found:
                break
        if found is None:
            continue
        case, samp = found
        rows.append({
            "file_id": h["file_id"], "file_name": h.get("file_name", ""),
            "file_size": h.get("file_size", 0),
            "patient_id": case["submitter_id"], "sample_id": samp["submitter_id"],
            "sample_type": samp["sample_type"], "group": SAMPLE_TYPES[samp["sample_type"]],
        })
    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError("No GDC files matched the query.")
    # one file per sample: keep the first file (sorted by file name) when a sample has several
    df = df.sort_values(["sample_id", "file_name"])
    dup = df["sample_id"].duplicated(keep="first")
    if dup.any():
        log.warning(f"{int(dup.sum())} additional file(s) for already-listed samples were dropped.")
    df = df[~dup].reset_index(drop=True)
    log.info(f"GDC returned {len(df)} samples: "
             + ", ".join(f"{k}={v}" for k, v in df['group'].value_counts().items()))
    return df


# ── Download and parse ─────────────────────────────────────────────────────────────────────────────

def download_file(file_id, dest, force=False, timeout=120, retries=4):
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0 and not force:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(GDC_DATA_URL + file_id, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
            tmp = dest.with_suffix(".part")
            tmp.write_bytes(data)
            tmp.replace(dest)
            return dest
        except Exception as exc:               # URLError, timeout, incomplete read, ...
            last = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"download of {file_id} failed after {retries} attempts: {last}")


def parse_star_counts(path):
    """Read one STAR counts file; returns a DataFrame indexed by gene_id (N_* summary rows removed)."""
    df = pd.read_csv(path, sep="\t", comment="#")
    df = df[~df["gene_id"].astype(str).str.startswith("N_")]
    return df.set_index("gene_id")


def build_matrices(sheet, files_dir, workers, force, protein_coding_only):
    """Download every file (thread pool) and assemble the counts / TPM matrices."""
    files_dir = Path(files_dir)
    t0 = time.time()
    paths = {}
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        futs = {pool.submit(download_file, r.file_id, files_dir / f"{r.file_id}.tsv", force): r
                for r in sheet.itertuples()}
        for i, fut in enumerate(as_completed(futs), 1):
            r = futs[fut]
            paths[r.sample_id] = fut.result()
            if i % 25 == 0 or i == len(futs):
                log.info(f"  downloaded {i}/{len(futs)} files ({time.time() - t0:.0f}s)")

    counts, tpm, annot = {}, {}, None
    for r in sheet.itertuples():
        t = parse_star_counts(paths[r.sample_id])
        if annot is None:
            annot = t[["gene_name", "gene_type"]]
        counts[r.sample_id] = t["unstranded"]
        tpm[r.sample_id] = t["tpm_unstranded"]
    counts = pd.DataFrame(counts)
    tpm = pd.DataFrame(tpm)

    if protein_coding_only:
        keep = annot.index[annot["gene_type"] == "protein_coding"]
        counts, tpm, annot = counts.loc[keep], tpm.loc[keep], annot.loc[keep]
        log.info(f"Restricted to {len(keep):,} protein-coding genes (use --all-genes to keep all).")

    # gene_id -> gene symbol; when several gene_ids share a symbol keep the one with the highest mean TPM
    symbols = annot["gene_name"].reindex(tpm.index)
    order = tpm.mean(axis=1).sort_values(ascending=False).index
    seen, selected = set(), []
    for g in order:
        s = symbols[g]
        if pd.notna(s) and s not in seen:                                # skip NaN / duplicate symbols
            seen.add(s)
            selected.append(g)
    counts = counts.loc[selected]
    tpm = tpm.loc[selected]
    counts.index = pd.Index(symbols[selected].to_numpy(), name="gene")
    tpm.index = pd.Index(symbols[selected].to_numpy(), name="gene")
    counts = counts.sort_index()
    tpm = tpm.sort_index()
    log.info(f"Matrices: {counts.shape[0]:,} genes x {counts.shape[1]:,} samples")
    return counts, tpm


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    default_out = Path(__file__).resolve().parents[1] / "data" / "raw" / "tcga_gdc"
    ap.add_argument("--out", default=str(default_out), help="output directory (default: data/raw/tcga_gdc)")
    ap.add_argument("--project", default="TCGA-LUAD")
    ap.add_argument("--workers", type=int, default=6, help="parallel downloads (default 6)")
    ap.add_argument("--paired-only", action="store_true",
                    help="keep only patients that have BOTH a primary tumor and a solid-tissue normal sample")
    ap.add_argument("--all-genes", action="store_true", help="keep non-protein-coding genes as well")
    ap.add_argument("--force", action="store_true", help="re-download cached files")
    ap.add_argument("--query-only", action="store_true", help="only list the files (no download)")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    sheet = query_files(args.project)
    if args.paired_only:
        both = set(sheet.loc[sheet.group == "tumor", "patient_id"]) & set(sheet.loc[sheet.group == "normal", "patient_id"])
        sheet = sheet[sheet["patient_id"].isin(both)].reset_index(drop=True)
        log.info(f"--paired-only: {len(both)} patients, {len(sheet)} samples kept.")
    sheet.to_csv(out / "gdc_files_manifest.csv", index=False)
    size_mb = sheet["file_size"].sum() / 1e6
    log.info(f"Total download size: about {size_mb:,.0f} MB")
    if args.query_only:
        return 0

    counts, tpm = build_matrices(sheet, out / "files", args.workers, args.force,
                                 protein_coding_only=not args.all_genes)
    counts.to_csv(out / "counts_matrix.csv")
    tpm.to_csv(out / "tpm_matrix.csv")
    sheet[["sample_id", "patient_id", "group", "sample_type", "file_id", "file_name"]].to_csv(
        out / "sample_sheet.csv", index=False)
    n_pairs = len(set(sheet.loc[sheet.group == "tumor", "patient_id"]) & set(sheet.loc[sheet.group == "normal", "patient_id"]))
    log.info(f"Wrote counts_matrix.csv, tpm_matrix.csv, sample_sheet.csv to {out}")
    log.info(f"Samples: {int((sheet.group == 'tumor').sum())} tumor, {int((sheet.group == 'normal').sum())} normal; "
             f"{n_pairs} patients with both.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
