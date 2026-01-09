import yaml


def validate_config(config_path="config.yaml"):
    with open(config_path) as f:
        cfg = yaml.safe_load(f)

    checks = []

    # Check multi-kmer pretrain settings
    checks.append(("multi_kmer_pretrain_enabled",
                   cfg.get("multi_kmer_pretrain_enabled") == True))
    checks.append(("pretrain_kmers",
                   cfg.get("pretrain_kmers") == [3, 4, 5, 6]))

    # Check early stopping settings
    checks.append(("mlm_early_stopping_enabled",
                   cfg.get("mlm_early_stopping_enabled") == True))
    checks.append(("mlm_early_stopping_patience",
                   cfg.get("mlm_early_stopping_patience", 0) > 0))

    # Check per-kmer overrides exist
    overrides = cfg.get("kmer_pretrain_overrides", {})
    for k in [3, 4, 5, 6]:
        checks.append((f"kmer_{k}_override_exists", k in overrides))

    # Check checkpoint paths
    mlm_ckpts = cfg.get("mlm_encoder_ckpt_by_k", {})
    for k in [3, 4, 5, 6]:
        checks.append((f"mlm_ckpt_k{k}_defined", k in mlm_ckpts))

    print("Multi-K-mer Config Validation")
    print("=" * 50)

    all_passed = True
    for name, passed in checks:
        status = "✓" if passed else "✗"
        print(f"{status} {name}")
        all_passed = all_passed and passed

    print("=" * 50)
    if all_passed:
        print("CONFIG VALIDATION PASSED")
    else:
        print("CONFIG VALIDATION FAILED")

    return all_passed


if __name__ == "__main__":
    validate_config()
