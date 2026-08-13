from pathlib import Path
import sys
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from measles_analysis.pipeline import run_pipeline


@pytest.fixture(scope="session")
def repo_root():
    return ROOT


@pytest.fixture(scope="session")
def generated_outputs(tmp_path_factory):
    out = tmp_path_factory.mktemp("full_pipeline")
    key = run_pipeline(ROOT / "data" / "derived", out, run_stochastic=True)
    return out, key
