"""Engineered features processing modules."""
import torch
import torch.nn as nn


class EngineeredFeatureMLP(nn.Module):
    """Separate MLP branch for processing engineered features.

    This gives engineered features (TNC, PseEIIP) their own
    representation space before fusion with sequence features.

    Includes optional pre-norm, gated hidden block, and residual projection.
    """

    def __init__(
        self,
        input_dim: int = 288,
        hidden_dim: int = 512,
        output_dim: int = 256,
        dropout: float = 0.3,
        use_pre_norm: bool = True,
        input_dropout: float = 0.05,
        use_gated: bool = True,
        use_residual: bool = True,
    ):
        super().__init__()
        self.use_pre_norm = use_pre_norm
        self.use_gated = use_gated
        self.use_residual = use_residual

        self.input_norm = nn.LayerNorm(input_dim) if use_pre_norm else None
        self.pre_norm_scale = (
            nn.Parameter(torch.tensor(0.0)) if use_pre_norm else None
        )
        self.input_dropout = nn.Dropout(input_dropout) if input_dropout > 0 else None

        self.mlp = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, output_dim),
            nn.LayerNorm(output_dim),
        )
        if use_gated:
            self.gate_proj = nn.Linear(hidden_dim, hidden_dim)
            self.gate_scale = nn.Parameter(torch.tensor(0.0))
        else:
            self.gate_proj = None
            self.gate_scale = None

        if use_residual:
            self.residual_proj = (
                nn.Identity()
                if input_dim == output_dim
                else nn.Linear(input_dim, output_dim)
            )
            self.residual_scale = nn.Parameter(torch.tensor(0.0))
        else:
            self.residual_proj = None
            self.residual_scale = None

        self.output_dim = output_dim

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: (B, input_dim) engineered features
        Returns:
            (B, output_dim) processed features
        """
        residual_in = x
        if self.input_norm is not None:
            normed = self.input_norm(x)
            scale = torch.tanh(self.pre_norm_scale)
            x = x + scale * (normed - x)

        if self.input_dropout is not None:
            x = self.input_dropout(x)

        # Process through MLP layers with optional gating
        x = self.mlp[0](x)
        x = self.mlp[1](x)
        x = self.mlp[2](x)
        x = self.mlp[3](x)

        if self.gate_proj is not None:
            gate = torch.sigmoid(self.gate_proj(x))
            scale = torch.tanh(self.gate_scale)
            x = x * (1.0 + scale * (gate - 0.5) * 2.0)

        x = self.mlp[4](x)
        x = self.mlp[5](x)
        x = self.mlp[6](x)
        x = self.mlp[7](x)
        x = self.mlp[8](x)
        x = self.mlp[9](x)

        if self.residual_proj is not None:
            residual = self.residual_proj(residual_in)
            scale = torch.tanh(self.residual_scale)
            x = x + scale * residual

        return x
