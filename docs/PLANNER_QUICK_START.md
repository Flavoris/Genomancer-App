# Multi-Step Cloning Planner - Quick Start

## Overview
The multi-step cloning planner is now fully functional after fixing the "unhashable type: Step" error. It can design complex cloning strategies automatically.

## Basic Usage

### 1. Minimal Test (Verify Installation)
```bash
python sim.py --plan-cloning tests/test_planner_spec.json --max-steps 2
```

Expected output:
```
✓ Plan found!
Total steps: 2
Score: 2.90
```

### 2. Run with Your Own Spec
Create a spec file (JSON or YAML):

```json
{
  "vector": {
    "name": "pVector",
    "fasta": "path/to/vector.fasta",
    "circular": true
  },
  "inserts": [
    {
      "name": "Insert1",
      "fasta": "path/to/insert.fasta",
      "features": [
        {
          "type": "CDS",
          "frame": 0,
          "start": 0,
          "end": 300,
          "label": "gene",
          "direction": "forward"
        }
      ]
    }
  ],
  "target": {
    "order": ["pVector", "Insert1"],
    "junctions": [
      {
        "left": "pVector",
        "right": "Insert1",
        "directional": true,
        "scar": "none",
        "keep_frame": false
      }
    ]
  },
  "constraints": {
    "avoid_internal_cuts": true,
    "dephosphorylate_backbone": true,
    "min_overhang": 4
  }
}
```

Then run:
```bash
python sim.py --plan-cloning my_spec.json --max-steps 3
```

## Command-Line Options

### Basic Options
- `--plan-cloning <file>`: Specification file (JSON or YAML)
- `--max-steps <N>`: Maximum cloning steps (default: 3)
- `--beam-width <N>`: Search beam width (default: 10)

### Enzyme Control
- `--avoid-enzymes "EcoRI,BamHI"`: Exclude specific enzymes
- `--allow-enzymes "BsaI,BsmBI"`: Restrict to specific enzymes
- `--prefer-typeIIS`: Prefer Golden Gate assemblies

### Output Control
- `--export-plan <dir>`: Export detailed protocols
- `--print-gels`: Include gel simulations

## Examples

### Classic Restriction Cloning
```bash
python sim.py --plan-cloning spec.json --max-steps 2
```

### Golden Gate Assembly
```bash
python sim.py --plan-cloning spec.json --prefer-typeIIS --max-steps 2
```

### Frame-Preserving Fusion
```bash
python sim.py --plan-cloning fusion.json --frame-check --max-steps 4
```

### Avoid Problematic Enzymes
```bash
python sim.py --plan-cloning spec.json --avoid-enzymes "EcoRI,NotI"
```

## Testing

### Run Regression Tests
```bash
# Specific regression test
python tests/test_planner_regression.py

# All planner tests
python -m pytest tests/test_planner*.py -v

# Full test suite
python -m pytest tests/ -v
```

### Direct Python Usage
```python
from planner import plan_from_spec
from sim import load_enzyme_database
import json

# Load enzyme database
enzyme_db = load_enzyme_database()

# Load spec
with open('spec.json', 'r') as f:
    spec = json.load(f)

# Run planner
plan = plan_from_spec(
    spec=spec,
    enzyme_db=enzyme_db,
    max_steps=3,
    options={'beam_width': 10}
)

print(f"Feasible: {plan.feasible}")
print(f"Steps: {len(plan.steps)}")
print(f"Score: {plan.score:.2f}")
```

## Troubleshooting

### "No feasible plan found"
- Increase `--max-steps`
- Relax constraints in spec
- Use `--allow-enzymes` to permit more enzymes
- Check that FASTA files exist

### "File not found"
- Use absolute paths or paths relative to working directory
- Verify FASTA files exist

### "Invalid specification"
- Check JSON/YAML syntax
- Ensure all required fields are present
- Validate FASTA file paths

## What Was Fixed

Previously, the planner crashed with:
```
TypeError: unhashable type: 'Step'
```

This has been fixed by implementing signature-based state hashing that doesn't require `Step` objects to be hashable. The planner now works correctly for all cloning scenarios.

## More Information

- Full documentation: `PLANNER_GUIDE.md`
- Fix details: `PLANNER_FIX_SUMMARY.md`
- Examples: `examples/` directory
- Tests: `tests/test_planner*.py`

## Support

If you encounter any issues:
1. Verify the planner works with the test spec: `python sim.py --plan-cloning tests/test_planner_spec.json --max-steps 2`
2. Run regression tests: `python tests/test_planner_regression.py`
3. Check that your spec file is valid JSON/YAML
4. Ensure enzyme database (`enzymes.json`) is present

The planner should now work reliably for all multi-step cloning scenarios!

