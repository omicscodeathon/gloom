"""
run_pipeline.py - Master Pipeline Runner
Executes all 18 steps of the LUAD ML Pipeline in sequence.
Usage: cd LUAD_ML_Pipeline && python run_pipeline.py
"""
import sys, time, logging
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
config.create_output_dirs()

logging.basicConfig(level=getattr(logging, config.LOG_LEVEL),
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(config.LOG_FILE), logging.StreamHandler(sys.stdout)])
log = logging.getLogger("run_pipeline")

STEPS = [
    ("Step 0  - Config",              "config",                          "create_output_dirs"),
    ("Step 1  - Data loading",        "step1_data_loading",              "run_data_loading"),
    ("Step 2  - Preprocessing",       "step2_preprocessing",             "run_preprocessing"),
    ("Step 3  - Harmonization",       "step3_harmonization",             "run_harmonization"),
    ("Step 4  - Differential expr",   "step4_differential_expression",   "run_differential_expression"),
    ("Step 5  - Expr features",       "step5_expression_features",       "run_expression_feature_construction"),
    ("Step 6  - Co-expr network",     "step6_coexpression_network",      "run_coexpression_network"),
    ("Step 7  - Network features",    "step7_network_features",          "run_network_feature_extraction"),
    ("Step 8  - Feature integration", "step8_feature_integration",       "run_feature_integration"),
    ("Step 9  - Label construction",  "step9_label_construction",        "run_label_construction"),
    ("Step 10 - Train/val split",     "step10_train_val_split",          "run_train_val_split"),
    ("Step 11 - Model training",      "step11_model_training",           "run_model_training"),
    ("Step 12 - Model evaluation",    "step12_model_evaluation",         "run_model_evaluation"),
    ("Step 13 - Feature importance",  "step13_feature_importance",       "run_feature_importance"),
    ("Step 14 - Gene ranking",        "step14_gene_ranking",             "run_gene_ranking"),
    ("Step 15 - Network annotation",  "step15_network_annotation",       "run_network_annotation"),
    ("Step 16 - Network export",      "step16_network_export",           "run_network_export"),
    ("Step 17 - Visualization",       "step17_interactive_visualization","run_interactive_visualization"),
    ("Step 18 - Final report",        "step18_final_report",             "run_final_report"),
]

def run_pipeline():
    log.info("="*70); log.info("LUAD ML PIPELINE - FULL RUN"); log.info("="*70)
    pipeline_start = time.time(); completed = []; failed_step = None
    for step_label, module_name, func_name in STEPS:
        log.info(f">>> {step_label}"); t0 = time.time()
        try:
            module = __import__(module_name)
            func   = getattr(module, func_name)
            func()
            elapsed = time.time() - t0
            log.info(f"<<< {step_label} DONE ({elapsed:.1f}s)")
            completed.append((step_label, elapsed, "OK"))
        except Exception as e:
            elapsed = time.time() - t0
            log.error(f"<<< {step_label} FAILED: {type(e).__name__}: {e}")
            completed.append((step_label, elapsed, f"FAILED: {e}"))
            failed_step = step_label; break
    total = time.time() - pipeline_start
    log.info("="*70); log.info("EXECUTION SUMMARY"); log.info("="*70)
    for label, elapsed, status in completed:
        icon = "OK" if status=="OK" else "FAIL"
        log.info(f"  [{icon}] {label:<45} {elapsed:>7.1f}s")
    log.info(f"  Total: {total:.1f}s ({total/60:.1f} min)")
    if failed_step: log.error(f"  STOPPED at: {failed_step}"); sys.exit(1)
    else: log.info("  Pipeline COMPLETE.")
    log.info("="*70)

if __name__ == "__main__":
    run_pipeline()
