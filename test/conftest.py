"""
pytest configuration: make the pure helper module ``gloom_utils`` importable.

The tests only import ``gloom_utils`` (no config, no file system side effects), and use
small synthetic data, so they run without any of the real LUAD data files.
"""
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
for candidate in (_REPO / "scripts", _REPO / "src" / "gloom" / "pipeline"):
    if (candidate / "gloom_utils.py").exists():
        sys.path.insert(0, str(candidate))
        break
