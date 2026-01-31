"""Pretraining components for Gene Whisperer MLM.

This package contains modular components for MLM pretraining:
- mlm_model.py: DNAMLM (Masked Language Model wrapper)

Note: Most pretraining functionality is still in pretrain_mlm.py
for backward compatibility. This module provides a cleaner import path
for key components.
"""

from .mlm_model import DNAMLM

__all__ = ["DNAMLM"]
