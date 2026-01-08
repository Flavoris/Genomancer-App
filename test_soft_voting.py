import torch

# Simulate outputs from 4 models (k=3,4,5,6)
# Model outputs are logits (before sigmoid)

# Scenario: Models disagree
model_logits = torch.tensor(
    [
        [[2.0]],  # k=3: sigmoid=0.88, confident positive
        [[1.0]],  # k=4: sigmoid=0.73, positive
        [[-0.5]],  # k=5: sigmoid=0.38, leaning negative
        [[0.5]],  # k=6: sigmoid=0.62, leaning positive
    ]
)

# Method 1: Average logits then sigmoid (WRONG)
avg_logits = model_logits.mean(dim=0)
wrong_prob = torch.sigmoid(avg_logits)
print(f"WRONG (avg logits then sigmoid): {wrong_prob.item():.4f}")

# Method 2: Sigmoid then average (CORRECT - soft voting)
probs = torch.sigmoid(model_logits)
correct_prob = probs.mean(dim=0)
print(f"CORRECT (sigmoid then avg): {correct_prob.item():.4f}")

# The difference matters!
print(f"Difference: {abs(wrong_prob.item() - correct_prob.item()):.4f}")

# With equal weights, soft voting gives: (0.88 + 0.73 + 0.38 + 0.62) / 4 = 0.6525
expected = 0.6525
if abs(correct_prob.item() - expected) < 0.01:
    print("SOFT VOTING TEST PASSED \u2713")
else:
    print(f"SOFT VOTING TEST FAILED \u2717 (expected {expected})")
