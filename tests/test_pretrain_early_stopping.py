from pathlib import Path
import sys

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from gene_whisperer.training.pretrain_mlm import _is_improvement


def test_is_improvement_respects_min_delta() -> None:
    assert _is_improvement(current_loss=1.0, best_loss=1.2, min_delta=0.0)
    assert _is_improvement(current_loss=1.0, best_loss=1.2, min_delta=0.1)
    assert not _is_improvement(current_loss=1.11, best_loss=1.2, min_delta=0.1)
    assert not _is_improvement(current_loss=1.2, best_loss=1.2, min_delta=0.0)
