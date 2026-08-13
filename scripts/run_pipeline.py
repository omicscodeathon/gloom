"""
run_pipeline.py — Master Pipeline Runner
=========================================
Executes the full LUAD ML Pipeline in sequence.

Usage:
  python run_pipeline.py                   # run all steps
  python run_pipeline.py --from 4          # resume from step 4 (skip 0-3)
  python run_pipeline.py --only 12         # run only step 12
  python run_pipeline.py --from 6 --to 11  # run steps 6 through 11

Step numbering:
  0   Config / output dirs
  1   Data loading
  1b  Batch correction (optional — requires USE_BATCH_CORRECTION=True)
  2   Preprocessing / QC
  3   Harmonization
  4   Differential expression
  5   Expression features
  6   Tumor co-expression network
  6b  Normal co-expression network
  7   Tumor network features
  7b  Differential network features
  8   Feature integration
  9   Label construction
  10  Train / validation split
  11  Model training
  11b PU Bagging (Mordelet-Vert)
  12  Model evaluation
  13  Feature importance
  14  Gene ranking
  15  Network annotation
  16  Network export
  17  Interactive visualization
  18  Final report
  19  KEGG enrichment
"""
import argparse, importlib, sys, time, logging
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()

logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE, encoding="utf-8"), logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger("run_pipeline")

# ── Step registry ──────────────────────────────────────────────────────────────
# Format: (step_key, display_label, module_name, function_name, optional)
#   optional=True  → warn on failure but continue
#   optional=False → halt pipeline on failure
STEPS = [
    ("00",  "Step 0   — Config",                "config",
             "create_output_dirs",                                          False),
    ("01",  "Step 1   — Data loading",          "step1_data_loading",
             "run_data_loading",                                            False),
    ("01b", "Step 1b  — Batch correction",      "step1b_batch_correction",
             "run_batch_correction",                                        True),
    ("02",  "Step 2   — Preprocessing",         "step2_preprocessing",
             "run_preprocessing",                                           False),
    ("03",  "Step 3   — Harmonization",         "step3_harmonization",
             "run_harmonization",                                           False),
    ("04",  "Step 4   — Differential expr",     "step4_differential_expression",
             "run_differential_expression",                                 False),
    ("05",  "Step 5   — Expr features",         "step5_expression_features",
             "run_expression_feature_construction",                         False),
    ("06",  "Step 6   — Tumor co-expr network", "step6_coexpression_network",
             "run_coexpression_network",                                    False),
    ("06b", "Step 6b  — Normal co-expr network","step6b_normal_network",
             "run_normal_coexpression_network",                             True),
    ("07",  "Step 7   — Network features",      "step7_network_features",
             "run_network_feature_extraction",                              False),
    ("07b", "Step 7b  — Diff. network features","step7b_differential_network_features",
             "run_differential_network_features",                           True),
    ("08",  "Step 8   — Feature integration",   "step8_feature_integration",
             "run_feature_integration",                                     False),
    ("09",  "Step 9   — Label construction",    "step9_label_construction",
             "run_label_construction",                                      False),
    ("10",  "Step 10  — Train/val split",       "step10_train_val_split",
             "run_train_val_split",                                         False),
    ("11",  "Step 11  — Model training",        "step11_model_training",
             "run_model_training",                                          False),
    ("11b", "Step 11b — PU Bagging",            "step11b_pu_bagging",
             "run_pu_bagging",                                              True),
    ("12",  "Step 12  — Model evaluation",      "step12_model_evaluation",
             "run_model_evaluation",                                        False),
    ("13",  "Step 13  — Feature importance",    "step13_feature_importance",
             "run_feature_importance",                                      False),
    ("14",  "Step 14  — Gene ranking",          "step14_gene_ranking",
             "run_gene_ranking",                                            False),
    ("15",  "Step 15  — Network annotation",    "step15_network_annotation",
             "run_network_annotation",                                      False),
    ("16",  "Step 16  — Network export",        "step16_network_export",
             "run_network_export",                                          False),
    ("17",  "Step 17  — Visualization",         "step17_interactive_visualization",
             "run_interactive_visualization",                               False),
    ("18",  "Step 18  — Final report",          "step18_final_report",
             "run_final_report",                                            False),
    ("19",  "Step 19  — KEGG enrichment",       "step19_kegg_enrichment",
             "run_kegg_enrichment",                                         True),
]


def _batch_correction_redirect_check():
    if not getattr(config, "USE_BATCH_CORRECTION", False):
        return
    corrected_tumor  = getattr(config, "BATCH_CORRECTED_TUMOR_FILE",  None)
    corrected_normal = getattr(config, "BATCH_CORRECTED_NORMAL_FILE", None)
    current_tumor    = config.TUMOR_EXPR_HARMONIZED
    if corrected_tumor and Path(str(corrected_tumor)).exists():
        if str(current_tumor) != str(corrected_tumor):
            log.warning(
                "\n" + "=" * 70 +
                "\n  BATCH CORRECTION REDIRECT REQUIRED\n"
                "  step1b wrote corrected files but config still points to the\n"
                "  un-corrected harmonized matrices. Update config.py:\n\n"
                f"    TUMOR_EXPR_HARMONIZED  = PROCESSED_DIR / "
                f'"{Path(str(corrected_tumor)).name}"\n'
                f"    NORMAL_EXPR_HARMONIZED = PROCESSED_DIR / "
                f'"{Path(str(corrected_normal)).name}"\n\n'
                "  Then re-run from step 4 for batch correction to take effect.\n"
                + "=" * 70
            )


def _filter_steps(steps, from_key=None, to_key=None, only_key=None, skip_optional=False):
    if only_key:
        filtered = [s for s in steps if s[0] == only_key]
        if skip_optional:
            filtered = [s for s in filtered if not s[4]]
        return filtered
    filtered = steps
    if from_key:
        from_idx = next((i for i, s in enumerate(steps) if s[0] == from_key), 0)
        filtered = filtered[from_idx:]
    if to_key:
        to_idx = next((i for i, s in enumerate(filtered) if filtered[i][0] == to_key), len(filtered) - 1)
        filtered = filtered[:to_idx + 1]
    if skip_optional:
        filtered = [s for s in filtered if not s[4]]
    return filtered


def run_pipeline(from_key=None, to_key=None, only_key=None, skip_optional=False):
    steps = _filter_steps(STEPS, from_key, to_key, only_key, skip_optional=skip_optional)
    if not steps:
        log.error("No steps selected to run.")
        sys.exit(2)

    log.info("=" * 70)
    log.info("LUAD ML PIPELINE — FULL RUN")
    log.info(f"Running {len(steps)} step(s)")
    if skip_optional:
        log.info("Optional steps are excluded for this run.")
    log.info("=" * 70)

    pipeline_start = time.time()
    completed      = []
    failed_step    = None

    for step_key, step_label, module_name, func_name, optional in steps:
        tag = "[optional]" if optional else ""
        log.info(f">>> {step_label} {tag}")
        t0 = time.time()
        try:
            module = importlib.import_module(module_name)
            importlib.reload(module)
            func = getattr(module, func_name)
            func()
            elapsed = time.time() - t0
            log.info(f"<<< {step_label} DONE ({elapsed:.1f}s)")
            completed.append((step_label, elapsed, "OK", optional))
        except Exception as e:
            elapsed = time.time() - t0
            if optional:
                log.warning(f"<<< {step_label} SKIPPED ({type(e).__name__}: {e})")
                completed.append((step_label, elapsed, f"SKIPPED: {e}", optional))
            else:
                log.error(f"<<< {step_label} FAILED: {type(e).__name__}: {e}")
                completed.append((step_label, elapsed, f"FAILED: {e}", optional))
                failed_step = step_label
                break

        if step_key == "01b":
            _batch_correction_redirect_check()

    total = time.time() - pipeline_start

    log.info("=" * 70)
    log.info("EXECUTION SUMMARY")
    log.info("=" * 70)
    for label, elapsed, status, optional in completed:
        if status == "OK":
            icon = " OK "
        elif status.startswith("SKIPPED"):
            icon = "SKIP"
        else:
            icon = "FAIL"
        opt_tag = " [opt]" if optional else "      "
        log.info(f"  [{icon}]{opt_tag} {label:<48} {elapsed:>7.1f}s")
    log.info(f"  Total: {total:.1f}s ({total/60:.1f} min)")

    if failed_step:
        log.error(f"  STOPPED at: {failed_step}")
        sys.exit(1)
    else:
        log.info("  Pipeline COMPLETE.")
    log.info("=" * 70)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="LUAD ML Pipeline runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Step keys:  0, 1, 1b, 2, 3, 4, 5, 6, 6b, 7, 7b, 8, 9, 10,\n"
            "            11, 11b, 12, 13, 14, 15, 16, 17, 18, 19\n\n"
            "Examples:\n"
            "  python run_pipeline.py                  # full run\n"
            "  python run_pipeline.py --skip-optional  # full run without optional steps\n"
            "  python run_pipeline.py --from 4         # resume from DE step\n"
            "  python run_pipeline.py --from 6 --to 8 # network to integration\n"
            "  python run_pipeline.py --only 12        # re-evaluate models only\n"
            "  python run_pipeline.py --only 11b       # PU bagging only\n"
        ),
    )
    parser.add_argument("--from",  dest="from_key",  default=None,
                        help="Start from this step key (inclusive)")
    parser.add_argument("--to",    dest="to_key",    default=None,
                        help="Stop after this step key (inclusive)")
    parser.add_argument("--only",  dest="only_key",  default=None,
                        help="Run only this single step key")
    parser.add_argument("--skip-optional", action="store_true",
                        help="Exclude all steps marked optional from the selected run")
    args = parser.parse_args()

    run_pipeline(
        from_key=args.from_key,
        to_key=args.to_key,
        only_key=args.only_key,
        skip_optional=args.skip_optional,
    )
