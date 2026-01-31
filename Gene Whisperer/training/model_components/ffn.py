"""Feed-forward network variants for Gene Whisperer models."""
import torch
import torch.nn as nn
import torch.nn.functional as F


class SwiGLU(nn.Module):
    """SwiGLU Feed-Forward Network.

    SwiGLU(x) = (Swish(W1·x) * W3·x) · W2

    Used in LLaMA, PaLM, and other modern LLMs.
    Better gradient flow than GELU.

    Args:
        in_features: Input dimension
        hidden_features: Hidden dimension (default: 8/3 * in_features, rounded to 64)
        out_features: Output dimension (default: in_features)
        bias: Whether to use bias in linear layers (default: False)
        dropout: Dropout probability (default: 0.0)
    """

    def __init__(
        self,
        in_features: int,
        hidden_features: int = None,
        out_features: int = None,
        bias: bool = False,
        dropout: float = 0.0,
    ):
        super().__init__()
        out_features = out_features or in_features
        hidden_features = hidden_features or int(in_features * 8 / 3)
        # Round to multiple of 64 for efficiency
        hidden_features = ((hidden_features + 63) // 64) * 64

        self.w1 = nn.Linear(in_features, hidden_features, bias=bias)
        self.w2 = nn.Linear(hidden_features, out_features, bias=bias)
        self.w3 = nn.Linear(in_features, hidden_features, bias=bias)
        self.dropout = nn.Dropout(dropout) if dropout > 0 else nn.Identity()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.dropout(self.w2(F.silu(self.w1(x)) * self.w3(x)))


class GLUFFN(nn.Module):
    """Gated Linear Unit Feed-Forward Network.

    Uses gating mechanism: output = Linear(GELU(gate) * up)
    where gate and up are separate linear projections.

    This is used in modern transformers like PaLM, LLaMA, and others
    for improved expressiveness over standard FFN.

    Note: GLU effectively doubles the parameters in the up-projection,
    so we use ff_dim * 2/3 for each branch to keep param count similar
    to standard FFN.
    """

    def __init__(
        self,
        d_model: int,
        ff_dim: int,
        dropout: float = 0.1,
        activation: str = "gelu",
    ):
        super().__init__()
        # To maintain similar param count, each branch gets 2/3 of ff_dim
        hidden_dim = int(ff_dim * 2 / 3)
        # Round to multiple of 8 for efficiency
        hidden_dim = ((hidden_dim + 7) // 8) * 8

        self.gate_proj = nn.Linear(d_model, hidden_dim, bias=False)
        self.up_proj = nn.Linear(d_model, hidden_dim, bias=False)
        self.down_proj = nn.Linear(hidden_dim, d_model, bias=False)

        self.dropout = nn.Dropout(dropout)

        if activation == "gelu":
            self.activation = nn.GELU()
        elif activation == "silu":
            self.activation = nn.SiLU()
        else:
            self.activation = nn.GELU()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: Input tensor of shape (batch, seq_len, d_model)
        Returns:
            Output tensor of shape (batch, seq_len, d_model)
        """
        gate = self.activation(self.gate_proj(x))
        up = self.up_proj(x)
        hidden = gate * up
        output = self.down_proj(hidden)
        return self.dropout(output)
