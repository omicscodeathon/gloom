"""
step20_candidate_evidence.py — External evidence for non-LCGene candidates (OPTIONAL, needs internet)
====================================================================================================
Absence from LCGene does not establish biological novelty (several highlighted genes such as
CDC20, MYBL2 or MARCO have published lung-cancer associations).  This step queries two public
resources for every non-LCGene gene in the top-N of the ranking (config.EVIDENCE_TOP_N, default 200)
and classifies it:

  already_LUAD             evidence for lung adenocarcinoma
  lung_cancer_unspecified  evidence for lung cancer / lung carcinoma without LUAD specificity
  other_cancer             evidence for cancer in general (other tumour types)
  indirect_only            some co-mention / weak association below the evidence thresholds
  no_association_found     BOTH sources queried successfully and returned nothing
  not_assessed             a source could not be reached and no positive evidence was found

ONLY ``no_association_found`` may be described as "potentially novel" — and only with respect to
these two sources and this search; it is a hypothesis for experimental follow-up, not a result.

Sources
  * Open Targets Platform GraphQL API  https://api.platform.opentargets.org/api/v4/graphql
      - gene symbol -> Ensembl id with the ``search`` query (entity = target)
      - ``target(ensemblId).associatedDiseases`` (indirect associations enabled, i.e. scores of
        descendant diseases propagate to their ancestors) restricted to
          lung adenocarcinoma      EFO_0000571
          lung carcinoma           EFO_0001071
          cancer (generic)         EFO_0000311 / MONDO_0004992
      The GraphQL schema evolves; several argument variants are tried in turn and the first that
      the server accepts is cached for the run.
  * Europe PMC REST  https://www.ebi.ac.uk/europepmc/webservices/rest/search
      hitCount of  "SYMBOL" AND "lung adenocarcinoma" / "lung cancer" / "cancer"
      (full-text co-mention counts, NOT curated associations; very short or common-word symbols
       can be ambiguous — flagged in ``short_symbol_flag``).

Robustness: timeouts, retries with exponential back-off, a local JSON cache
(results/candidate_evidence_cache.json; delete it to refresh) and a connectivity probe.  The step
never raises into the pipeline: offline or on any error it logs a warning and returns {}.

Output: results/candidate_evidence.csv
"""
import json
import logging
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
from gloom_utils import classify_candidate_evidence

config.create_output_dirs()
logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE, encoding="utf-8"), logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

OT_URL   = "https://api.platform.opentargets.org/api/v4/graphql"
EPMC_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
USER_AGENT = "gloom-candidate-evidence/0.2 (research; https://github.com/omicscodeathon/gloom)"

EFO_LUAD = "EFO_0000571"
EFO_LUNG_CARCINOMA = "EFO_0001071"
EFO_CANCER_IDS = ("EFO_0000311", "MONDO_0004992")
DISEASE_FILTER = [EFO_LUAD, EFO_LUNG_CARCINOMA, *EFO_CANCER_IDS]

# ── GraphQL documents (variants for different Open Targets schema versions) ────────────────────────
Q_SEARCH = """
query SearchTarget($q: String!) {
  search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 10}) {
    hits { id name entity }
  }
}
"""

Q_ASSOC_VARIANTS = [
    # current API: Bs filter + enableIndirect + page
    """
    query Assoc($id: String!, $bs: [String!]) {
      target(ensemblId: $id) {
        id approvedSymbol
        associatedDiseases(Bs: $bs, enableIndirect: true, page: {index: 0, size: 50}) {
          count
          rows { disease { id name } score }
        }
      }
    }
    """,
    # same without paging arguments
    """
    query Assoc($id: String!, $bs: [String!]) {
      target(ensemblId: $id) {
        id approvedSymbol
        associatedDiseases(Bs: $bs, enableIndirect: true) {
          count
          rows { disease { id name } score }
        }
      }
    }
    """,
    # older API: index/size arguments
    """
    query Assoc($id: String!, $bs: [String!]) {
      target(ensemblId: $id) {
        id approvedSymbol
        associatedDiseases(Bs: $bs, enableIndirect: true, index: 0, size: 50) {
          count
          rows { disease { id name } score }
        }
      }
    }
    """,
]

_cache_lock = threading.Lock()
_variant_lock = threading.Lock()
_working_variant = {"idx": None}


# ── HTTP helper ───────────────────────────────────────────────────────────────────────────────────

def _http_json(url, payload=None, timeout=20, retries=3, backoff=1.6):
    """GET (payload=None) or POST JSON; retries on network errors / 5xx / 429 with back-off.
    A 400 response carrying JSON (GraphQL validation errors) is returned, not retried."""
    last = None
    for attempt in range(max(1, retries)):
        try:
            headers = {"Accept": "application/json", "User-Agent": USER_AGENT}
            data = None
            if payload is not None:
                data = json.dumps(payload).encode("utf-8")
                headers["Content-Type"] = "application/json"
            req = urllib.request.Request(url, data=data, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            last = exc
            if exc.code == 400:
                try:
                    return json.loads(exc.read().decode("utf-8"))
                except Exception:
                    break
            if exc.code in (401, 403, 404):
                break
        except Exception as exc:                      # URLError, timeout, JSON error, ...
            last = exc
        time.sleep(backoff ** attempt)
    raise RuntimeError(f"request failed after {retries} attempt(s): {last}")


def _graphql(query, variables, timeout, retries):
    return _http_json(OT_URL, {"query": query, "variables": variables}, timeout, retries)


# ── Open Targets ──────────────────────────────────────────────────────────────────────────────────

def _resolve_ensembl(symbol, timeout, retries):
    """Gene symbol -> Ensembl gene id through the search query (exact symbol match preferred)."""
    resp = _graphql(Q_SEARCH, {"q": symbol}, timeout, retries)
    if resp.get("errors"):
        raise RuntimeError(f"search: {resp['errors'][0].get('message', resp['errors'])}")
    hits = ((resp.get("data") or {}).get("search") or {}).get("hits") or []
    hits = [h for h in hits if h.get("entity", "target") == "target" and str(h.get("id", "")).startswith("ENSG")]
    for h in hits:
        if str(h.get("name", "")).upper() == symbol.upper():
            return h["id"]
    return None          # no exact symbol match: do not guess a different gene


def _ot_scores(ensembl_id, timeout, retries):
    """Association scores for LUAD, lung carcinoma and generic cancer (max over cancer ids)."""
    variables = {"id": ensembl_id, "bs": DISEASE_FILTER}
    with _variant_lock:
        order = list(range(len(Q_ASSOC_VARIANTS)))
        if _working_variant["idx"] is not None:
            order.remove(_working_variant["idx"])
            order.insert(0, _working_variant["idx"])
    last_err = None
    for idx in order:
        resp = _graphql(Q_ASSOC_VARIANTS[idx], variables, timeout, retries)
        if resp.get("errors") or not (resp.get("data") or {}).get("target"):
            last_err = (resp.get("errors") or [{"message": "empty target"}])[0].get("message")
            continue
        with _variant_lock:
            _working_variant["idx"] = idx
        rows = (resp["data"]["target"].get("associatedDiseases") or {}).get("rows") or []
        by_id = {r["disease"]["id"]: float(r["score"]) for r in rows if r.get("disease")}
        return {
            "luad": by_id.get(EFO_LUAD, 0.0),
            "lung_carcinoma": by_id.get(EFO_LUNG_CARCINOMA, 0.0),
            "cancer": max([by_id.get(i, 0.0) for i in EFO_CANCER_IDS] or [0.0]),
        }
    raise RuntimeError(f"associatedDiseases query rejected by the server: {last_err}")


def query_open_targets(symbol, timeout, retries):
    try:
        ensembl = _resolve_ensembl(symbol, timeout, retries)
        if ensembl is None:
            return {"status": "no_target_in_open_targets", "ensembl_id": None,
                    "luad": 0.0, "lung_carcinoma": 0.0, "cancer": 0.0}
        s = _ot_scores(ensembl, timeout, retries)
        return {"status": "ok", "ensembl_id": ensembl, **s}
    except Exception as exc:
        return {"status": f"error: {exc}"[:200], "ensembl_id": None,
                "luad": None, "lung_carcinoma": None, "cancer": None}


# ── Europe PMC ────────────────────────────────────────────────────────────────────────────────────

def _epmc_hits(symbol, phrase, timeout, retries):
    q = f'"{symbol}" AND "{phrase}"'
    url = (f"{EPMC_URL}?query={urllib.parse.quote(q)}&format=json&resultType=lite&pageSize=1")
    resp = _http_json(url, None, timeout, retries)
    return int(resp.get("hitCount", 0))


def query_europe_pmc(symbol, timeout, retries):
    out = {"status": "ok"}
    for key, phrase in (("luad", "lung adenocarcinoma"), ("lung_cancer", "lung cancer"), ("cancer", "cancer")):
        try:
            out[key] = _epmc_hits(symbol, phrase, timeout, retries)
        except Exception as exc:
            out[key] = None
            out["status"] = f"error: {exc}"[:200]
    return out


# ── Cache ─────────────────────────────────────────────────────────────────────────────────────────

def _load_cache(path):
    if Path(path).exists():
        try:
            return json.loads(Path(path).read_text(encoding="utf-8"))
        except Exception:
            log.warning(f"  Cache {path} unreadable — starting a new cache.")
    return {}


def _save_cache(path, cache):
    with _cache_lock:
        tmp = Path(str(path) + ".tmp")
        tmp.write_text(json.dumps(cache, indent=1), encoding="utf-8")
        tmp.replace(path)


def _cache_ok(rec):
    """A cached record is reused only if both sources answered."""
    return bool(rec) and rec.get("ot", {}).get("status", "").split(":")[0] in ("ok", "no_target_in_open_targets") \
        and rec.get("epmc", {}).get("status") == "ok"


def _fetch_gene(symbol, cache, timeout, retries):
    rec = cache.get(symbol)
    if _cache_ok(rec):
        return symbol, rec
    rec = {
        "ot": query_open_targets(symbol, timeout, retries),
        "epmc": query_europe_pmc(symbol, timeout, retries),
        "fetched": str(date.today()),
    }
    time.sleep(0.05)
    return symbol, rec


def _online(timeout=8):
    """Connectivity probe against Europe PMC (one cheap request, no retries)."""
    try:
        _http_json(f"{EPMC_URL}?query=%22TP53%22&format=json&resultType=lite&pageSize=1",
                   None, timeout, retries=1)
        return True
    except Exception as exc:
        log.warning(f"  Connectivity probe failed ({exc}).")
        return False


# ── Main ──────────────────────────────────────────────────────────────────────────────────────────

def run_candidate_evidence(force: bool = False) -> dict:
    """Optional step: returns {} (and never raises) when skipped, offline or on error."""
    log.info("=" * 60)
    log.info("STEP 20 — EXTERNAL EVIDENCE FOR NON-LCGENE CANDIDATES (optional)")
    log.info("=" * 60)
    try:
        if not force and not getattr(config, "USE_CANDIDATE_EVIDENCE", True):
            log.info("  USE_CANDIDATE_EVIDENCE = False — skipped.")
            return {}
        if not Path(config.GENE_RANKINGS_FILE).exists():
            log.warning(f"  {config.GENE_RANKINGS_FILE} not found — run step14 first. Skipped.")
            return {}
        if not _online():
            log.warning("  No internet access to Open Targets / Europe PMC — step20 skipped "
                        "(the main results are unaffected). Re-run this step with internet access.")
            return {}

        top_n = int(getattr(config, "EVIDENCE_TOP_N", 200))
        timeout = float(getattr(config, "EVIDENCE_TIMEOUT_S", 20))
        retries = int(getattr(config, "EVIDENCE_RETRIES", 3))
        workers = int(getattr(config, "EVIDENCE_WORKERS", 4))
        min_hits = int(getattr(config, "EVIDENCE_MIN_EPMC_HITS", 3))
        min_ot = float(getattr(config, "EVIDENCE_MIN_OT_SCORE", 0.05))

        ranking = pd.read_csv(config.GENE_RANKINGS_FILE, index_col=0)
        cand = ranking[~ranking["is_lcgene_gene"].astype(bool)].sort_values("rank").head(top_n)
        symbols = [str(g).strip().upper() for g in cand.index]
        log.info(f"  Checking the top {len(symbols)} non-LCGene genes (workers={workers}) …")

        cache_path = Path(config.RESULTS_DIR) / "candidate_evidence_cache.json"
        cache = _load_cache(cache_path)
        results = {}
        t0 = time.time()
        with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
            futures = {pool.submit(_fetch_gene, s, cache, timeout, retries): s for s in symbols}
            for i, fut in enumerate(as_completed(futures), 1):
                try:
                    sym, rec = fut.result()
                except Exception as exc:                 # defensive: _fetch_gene already catches
                    sym = futures[fut]
                    rec = {"ot": {"status": f"error: {exc}", "luad": None, "lung_carcinoma": None, "cancer": None},
                           "epmc": {"status": f"error: {exc}", "luad": None, "lung_cancer": None, "cancer": None},
                           "fetched": str(date.today())}
                results[sym] = rec
                cache[sym] = rec
                if i % 25 == 0 or i == len(futures):
                    _save_cache(cache_path, cache)
                    log.info(f"  {i}/{len(futures)} genes queried ({time.time() - t0:.0f}s)")

        rows = []
        for gene, sym in zip(cand.index, symbols):
            rec, row = results[sym], cand.loc[gene]
            ot, ep = rec["ot"], rec["epmc"]
            cls = classify_candidate_evidence(
                ot.get("luad"), ot.get("lung_carcinoma"), ot.get("cancer"),
                ep.get("luad"), ep.get("lung_cancer"), ep.get("cancer"),
                min_ot_score=min_ot, min_hits=min_hits)
            rows.append({
                "gene": gene,
                "rank": int(row["rank"]),
                "predicted_prob": float(row["predicted_prob"]),
                "log2fc": float(row.get("log2fc", np.nan)),
                "direction": row.get("direction", ""),
                "ensembl_id": ot.get("ensembl_id"),
                "ot_status": ot.get("status"),
                "ot_luad_score": ot.get("luad"),
                "ot_lung_carcinoma_score": ot.get("lung_carcinoma"),
                "ot_cancer_score": ot.get("cancer"),
                "epmc_status": ep.get("status"),
                "epmc_luad_hits": ep.get("luad"),
                "epmc_lung_cancer_hits": ep.get("lung_cancer"),
                "epmc_cancer_hits": ep.get("cancer"),
                "short_symbol_flag": len(sym) <= 3,
                "evidence_class": cls,
                "potentially_novel": cls == "no_association_found",
            })
        out = pd.DataFrame(rows)
        out_path = Path(config.RESULTS_DIR) / "candidate_evidence.csv"
        out.to_csv(out_path, index=False)
        log.info("\n" + out["evidence_class"].value_counts().to_string())
        log.info(f"  Saved -> {out_path}")
        log.info("  Reminder: only 'no_association_found' may be called 'potentially novel' "
                 "(Open Targets + Europe PMC only; co-mention counts are not curated associations).")
        log.info("STEP 20 COMPLETE")
        return {"evidence": out}
    except Exception as exc:                              # never break the main pipeline
        log.warning(f"  Step 20 aborted ({type(exc).__name__}: {exc}); main results unaffected.")
        return {}


if __name__ == "__main__":
    r = run_candidate_evidence(force=True)
    if r:
        print(r["evidence"]["evidence_class"].value_counts())
