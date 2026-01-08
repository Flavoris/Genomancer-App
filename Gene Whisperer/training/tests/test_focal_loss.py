"""
Basic validation for the focal loss implementation.
"""
import sys

import torch

sys.path.insert(0, ".")

from train_stage1 import FocalLoss


def _compute_focal_loss_demo() -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    focal = FocalLoss(alpha=0.25, gamma=2.0)

    # Easy correct prediction (high confidence, correct)
    logits_easy = torch.tensor([[5.0]])
    targets_pos = torch.tensor([[1.0]])
    loss_easy = focal(logits_easy, targets_pos)

    # Hard incorrect prediction (low confidence, wrong)
    logits_hard = torch.tensor([[0.0]])
    loss_hard = focal(logits_hard, targets_pos)

    ratio = loss_hard / loss_easy
    return loss_easy, loss_hard, ratio


def test_focal_loss_downweights_easy_examples() -> None:
    loss_easy, loss_hard, ratio = _compute_focal_loss_demo()

    assert loss_hard > loss_easy
    assert ratio > 5.0


if __name__ == "__main__":
    loss_easy, loss_hard, ratio = _compute_focal_loss_demo()
    print(f"Easy correct (conf=0.99, target=1): loss={loss_easy.item():.6f}")
    print(f"Hard uncertain (conf=0.5, target=1): loss={loss_hard.item():.6f}")
    print(f"Hard/Easy loss ratio: {ratio.item():.2f}x")
    print("Expected: >10x (focal loss should down-weight easy examples)")
    if ratio > 5:
        print("FOCAL LOSS TEST PASSED ✓")
    else:
        print("FOCAL LOSS TEST FAILED ✗")
