from pathlib import Path
import sys

import numpy as np

ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from gene_whisperer.features.engineered import FeatureExtractor


def test_feature_extractor_dimensions() -> None:
    positives = ["ACGTA", "ACGTT", "ACGTC"]
    negatives = ["TTGCA", "TTGCT", "TTGCC"]

    extractor = FeatureExtractor(use_pstnp=True)
    extractor.fit_pstnp(positives, negatives)

    features = extractor.transform("ACGTA")
    assert features.shape[0] == extractor.feature_dim()
    assert np.isfinite(features).all()


def test_feature_extractor_gc_content() -> None:
    extractor = FeatureExtractor(use_pstnp=False, use_cksnap=False, use_tnc=False, use_pseeiip=False)
    features = extractor.transform("GGCC")
    assert features.shape[0] == extractor.feature_dim()
    gc, at = features.tolist()
    assert gc == 1.0
    assert at == 0.0
