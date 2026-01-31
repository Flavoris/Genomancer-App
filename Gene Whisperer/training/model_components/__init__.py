"""Gene Whisperer Model Components.

This package contains modular components for the Gene Whisperer model:
- norm.py: Normalization layers (RMSNorm)
- ffn.py: Feed-forward network variants (SwiGLU, GLUFFN)
- layers.py: Transformer layers (PreNormTransformerLayer, StochasticDepth, etc.)
- encoder.py: DNA encoder (DNAEncoder, AttentionPooling)
- features.py: Engineered features processing (EngineeredFeatureMLP)
- classifier.py: Stage 1 classifier models
"""

# Normalization
from .norm import RMSNorm

# Feed-forward networks
from .ffn import SwiGLU, GLUFFN

# Core layers
from .layers import (
    StochasticDepth,
    RelativePositionBias,
    PreNormTransformerLayer,
)

# Encoder and attention
from .encoder import (
    DNAEncoder,
    AttentionPooling,
    _load_checkpoint_file,
)

# Engineered features
from .features import EngineeredFeatureMLP

# Stage 1 classifiers
from .classifier import (
    GeneWhispererStage1,
    GeneWhispererTransformerOnly,
    GeneWhispererFeaturesOnly,
    GeneWhispererCombined,
    _mean_pool_sequence,
)

__all__ = [
    # Normalization
    "RMSNorm",
    # FFN
    "SwiGLU",
    "GLUFFN",
    # Layers
    "StochasticDepth",
    "RelativePositionBias",
    "PreNormTransformerLayer",
    # Encoder
    "DNAEncoder",
    "AttentionPooling",
    "_load_checkpoint_file",
    # Features
    "EngineeredFeatureMLP",
    # Classifiers
    "GeneWhispererStage1",
    "GeneWhispererTransformerOnly",
    "GeneWhispererFeaturesOnly",
    "GeneWhispererCombined",
    "_mean_pool_sequence",
]
