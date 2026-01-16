from __future__ import annotations

from ablation_variants import PIPELINE_VARIANTS, STAGE1_VARIANTS, STAGE2_VARIANTS


def test_stage_variants_include_all_engineered_groups() -> None:
    stage1_names = {name for name, _ in STAGE1_VARIANTS}
    stage2_names = {name for name, _ in STAGE2_VARIANTS}

    for name in ("no_tnc", "no_pseeiip", "no_cksnap", "no_pstnp"):
        assert name in stage1_names
        assert name in stage2_names


def test_pipeline_variants_include_baseline() -> None:
    pipeline_names = [name for name, _ in PIPELINE_VARIANTS]
    assert pipeline_names[0] == "baseline"
