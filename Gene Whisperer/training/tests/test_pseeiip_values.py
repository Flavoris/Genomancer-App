"""
Unit tests for PseEIIP feature calculation.
"""
import sys

sys.path.insert(0, ".")


def test_pseeiip_single_trinucleotide_uses_eiip_sum():
    """AAA should map to the EIIP sum of three A bases with frequency 1."""
    from dataset import BASE_TO_IDX, EIIP_VALUES, compute_pseeiip

    features = compute_pseeiip("AAA")

    aaa_idx = BASE_TO_IDX["A"] * 16 + BASE_TO_IDX["A"] * 4 + BASE_TO_IDX["A"]
    expected_value = EIIP_VALUES["A"] * 3.0
    actual_value = features[aaa_idx].item()

    assert features.shape[0] == 64
    assert abs(actual_value - expected_value) < 1e-6
    assert (features != 0).sum().item() == 1
    assert abs(features.sum().item() - expected_value) < 1e-6
    assert features.min().item() >= 0.0
