# Cloning Planner Guide

## Overview

The Genomancer cloning planner is an intelligent tool that automatically designs multi-step molecular cloning strategies. Given a target construct specification, it searches for an optimal sequence of digests, ligations, and assemblies to achieve your desired final construct.

## Features

- **Automated Strategy Design**: Searches through possible enzyme combinations and ligation orders
- **Multiple Cloning Methods**: Supports classic restriction cloning, Golden Gate assembly, and PCR-based approaches
- **Constraint-Aware**: Respects biological constraints like reading frame preservation, internal site avoidance, and directionality
- **Intelligent Scoring**: Ranks plans based on complexity, directional safety, and laboratory practicality
- **Export Protocols**: Generates detailed lab protocols and intermediate construct files

## Quick Start

### 1. Create a Specification File

Create a JSON or YAML file describing your cloning project:

```json
{
  "vector": {
    "name": "pBackbone",
    "fasta": "plasmid.fasta",
    "circular": true
  },
  "inserts": [
    {
      "name": "GeneA",
      "fasta": "geneA.fasta",
      "features": [
        {"type": "CDS", "frame": 0, "start": 0, "end": 300}
      ]
    }
  ],
  "target": {
    "order": ["pBackbone", "GeneA"],
    "junctions": [
      {
        "left": "pBackbone",
        "right": "GeneA",
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

### 2. Run the Planner

```bash
python sim.py --plan-cloning my_cloning.json --max-steps 3
```

### 3. View the Results

The planner will output:
- A summary of the optimal strategy found
- Step-by-step protocol details
- Predicted fragment sizes
- Score and feasibility assessment

## Specification Format

### Required Fields

#### `vector`
The plasmid backbone for your cloning.

```json
"vector": {
  "name": "pUC19",           // Construct name
  "fasta": "path/to/file",   // Path to FASTA file
  "circular": true           // Is it circular?
}
```

#### `target`
The desired final assembly.

```json
"target": {
  "order": ["vector", "insert1", "insert2"],  // Assembly order
  "junctions": [
    {
      "left": "vector",
      "right": "insert1",
      "directional": true,    // Require directional ligation?
      "scar": "none",         // Expected scar sequence
      "keep_frame": false     // Preserve reading frame?
    }
  ]
}
```

### Optional Fields

#### `inserts`
List of DNA fragments to insert.

```json
"inserts": [
  {
    "name": "GeneA",
    "fasta": "geneA.fasta",
    "features": [
      {
        "type": "CDS",
        "frame": 0,
        "start": 0,
        "end": 720,
        "label": "gfp",
        "direction": "forward"
      }
    ]
  }
]
```

#### `constraints`
Planning constraints.

```json
"constraints": {
  "avoid_internal_cuts": true,       // Don't cut inside CDS features
  "dephosphorylate_backbone": true,  // DeP the vector backbone
  "min_overhang": 4                  // Minimum overhang length (bp)
}
```

## Command-Line Options

### Basic Options

- `--plan-cloning <file>`: Specify the spec file (JSON or YAML)
- `--max-steps <N>`: Maximum number of cloning steps (default: 3)
- `--export-plan <dir>`: Export detailed protocols to directory
- `--print-gels`: Include simulated gel images for each step

### Enzyme Constraints

- `--avoid-enzymes "EcoRI,NotI"`: Forbid specific enzymes
- `--allow-enzymes "BsaI,SpeI,XbaI"`: Restrict to specific enzymes (whitelist)
- `--prefer-typeIIS`: Bias toward Golden Gate assemblies

### Advanced Options

- `--frame-check`: Enforce reading-frame preservation at junctions
- `--require-directional`: Only consider directional ligation strategies

## Planning Strategies

### Classic Restriction Cloning

The planner can design traditional two-enzyme directional cloning:

```bash
python sim.py --plan-cloning spec.json --max-steps 2
```

**Example output:**
```
Step 1 — Digest: pBackbone with EcoRI + BamHI (deP: YES)
  Fragments: [backbone 4123 bp, insert_window 1276 bp]
  
Step 2 — Ligate: backbone + GeneA
  Junction: directional ✔; scar: none
```

### Golden Gate Assembly

For multi-part assemblies, enable Type IIS preference:

```bash
python sim.py --plan-cloning spec.json --prefer-typeIIS --max-steps 2
```

**Example output:**
```
Step 1 — Golden Gate (BsaI): assemble 3 parts
  Overhangs: AATG→AGGT→GCTT
  
Final — Construct: pBackbone+GeneA+GeneB, 6210 bp (circular)
```

### Frame-Preserving Fusions

For protein fusions, enable frame checking:

```bash
python sim.py --plan-cloning fusion_spec.json --frame-check --max-steps 4
```

## Scoring System

Plans are scored based on multiple factors (lower scores are better):

| Factor | Weight | Description |
|--------|--------|-------------|
| Steps | 1.0 | Number of cloning steps |
| Unique enzymes | 0.5 | Different enzymes used |
| Buffer switches | 0.3 | Estimated buffer changes |
| Non-directional | 1.0 | Penalty if directionality lost |
| Internal cuts | 2.0 | Cutting inside protected features |
| Scar length | 0.1 | Length of junction scars |
| Type IIS bonus | -0.4 | Reward for Golden Gate |
| Enzyme reuse | -0.3 | Reward for using same enzyme |

## Examples

### Example 1: Simple Gene Insertion

```json
{
  "vector": {"name": "pUC19", "fasta": "puc19.fasta", "circular": true},
  "inserts": [{"name": "GFP", "fasta": "gfp.fasta"}],
  "target": {
    "order": ["pUC19", "GFP"],
    "junctions": [{"left": "pUC19", "right": "GFP", "directional": true}]
  },
  "constraints": {"avoid_internal_cuts": true}
}
```

**Command:**
```bash
python sim.py --plan-cloning insert_gfp.json --max-steps 2
```

### Example 2: Three-Part Golden Gate

```json
{
  "vector": {"name": "pDest", "fasta": "dest.fasta", "circular": true},
  "inserts": [
    {"name": "Promoter", "fasta": "promoter.fasta"},
    {"name": "CDS", "fasta": "cds.fasta"},
    {"name": "Terminator", "fasta": "term.fasta"}
  ],
  "target": {
    "order": ["pDest", "Promoter", "CDS", "Terminator"],
    "junctions": [
      {"left": "pDest", "right": "Promoter", "scar": "AATG"},
      {"left": "Promoter", "right": "CDS", "scar": "AGGT"},
      {"left": "CDS", "right": "Terminator", "scar": "GCTT"},
      {"left": "Terminator", "right": "pDest", "scar": "TCGA"}
    ]
  }
}
```

**Command:**
```bash
python sim.py --plan-cloning gg_assembly.json --prefer-typeIIS --export-plan gg_protocol/
```

### Example 3: Fusion Protein with Frame Check

```json
{
  "vector": {"name": "pExp", "fasta": "expression.fasta", "circular": true},
  "inserts": [
    {
      "name": "GFP",
      "fasta": "gfp.fasta",
      "features": [{"type": "CDS", "frame": 0, "start": 0, "end": 720}]
    },
    {
      "name": "Linker_mCherry",
      "fasta": "linker_mcherry.fasta",
      "features": [{"type": "CDS", "frame": 0, "start": 15, "end": 735}]
    }
  ],
  "target": {
    "order": ["pExp", "GFP", "Linker_mCherry"],
    "junctions": [
      {"left": "pExp", "right": "GFP", "keep_frame": true},
      {"left": "GFP", "right": "Linker_mCherry", "keep_frame": true}
    ]
  }
}
```

**Command:**
```bash
python sim.py --plan-cloning fusion.json --frame-check --max-steps 4
```

## Export Formats

When using `--export-plan`, the planner generates:

- `plan.json`: Machine-readable plan specification
- `step_01_digest.gb`: GenBank file for step 1
- `step_01_fragments.csv`: Fragment data for step 1
- `step_02_ligate.gb`: GenBank file for step 2
- `final.gb`: Final construct GenBank file

## Troubleshooting

### "No feasible plan found"

**Possible causes:**
- No suitable enzyme combinations exist
- Constraints are too restrictive
- Max steps is too low

**Solutions:**
- Increase `--max-steps`
- Relax constraints in spec
- Check that enzymes cut at appropriate sites
- Use `--allow-enzymes` to specify known working enzymes

### "Internal cuts detected"

**Solution:**
Use `--avoid-enzymes` to exclude problematic enzymes:
```bash
python sim.py --plan-cloning spec.json --avoid-enzymes "EcoRI,BamHI"
```

### Frame shift in fusion

**Solution:**
Ensure `keep_frame: true` in junctions and use `--frame-check`:
```bash
python sim.py --plan-cloning spec.json --frame-check
```

## Advanced Topics

### Custom Enzyme Database

The planner uses `enzymes.json` from the Genomancer directory. You can add custom enzymes:

```json
{
  "name": "MyEnzyme",
  "site": "GATATC",
  "cut_index": 3,
  "overhang_type": "Blunt"
}
```

### Beam Search Tuning

The planner uses beam search with default width of 10. For complex problems, you can modify this in the code:

```python
planner_options['beam_width'] = 20  # Wider search
```

### Heuristic Customization

Scoring weights can be adjusted in `planner.py`:

```python
def score_plan(plan: Plan, options: Dict) -> float:
    score = 0.0
    score += 1.0 * len(plan.steps)        # Adjust weight here
    score += 0.5 * len(enzymes_used)      # Adjust weight here
    # ...
```

## API Usage

You can also use the planner programmatically:

```python
from planner import plan_from_spec
from planner_utils import load_json_or_yaml

# Load specification
spec = load_json_or_yaml('my_cloning.json')

# Load enzyme database
from sim import load_enzyme_database
enzyme_db = load_enzyme_database()

# Run planner
plan = plan_from_spec(
    spec=spec,
    enzyme_db=enzyme_db,
    max_steps=3,
    options={'prefer_typeIIS': True}
)

# Check results
if plan.feasible:
    print(f"Found plan with {len(plan.steps)} steps")
    print(f"Score: {plan.score:.2f}")
else:
    print(f"No plan found: {plan.reason}")
```

## Citation

If you use the Genomancer cloning planner in your research, please cite:

```
Genomancer: An intelligent molecular cloning planner
[Your citation information here]
```

## Contributing

Contributions welcome! Key areas for improvement:

- Additional action types (Gibson assembly, USER cloning, etc.)
- Better heuristics for plan quality
- Support for more complex constraints
- Improved Golden Gate overhang design

## License

[Your license information here]

