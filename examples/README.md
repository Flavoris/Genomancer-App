# Cloning Planner Examples

This directory contains example specification files for the multi-step cloning planner.

## Specification Format

Specifications can be written in JSON or YAML format. Each spec contains:

- **vector**: The plasmid backbone
  - `name`: Construct name
  - `fasta`: Path to FASTA file
  - `circular`: Boolean indicating if circular

- **inserts**: List of DNA fragments to insert
  - `name`: Insert name
  - `fasta`: Path to FASTA file
  - `features`: Optional list of features (CDS, promoters, etc.)

- **target**: Desired final assembly
  - `order`: List of construct names in order
  - `junctions`: List of junction specifications
    - `left`: Name of left construct
    - `right`: Name of right construct
    - `directional`: Require directional ligation
    - `scar`: Expected scar sequence (or "none")
    - `keep_frame`: Preserve reading frame

- **constraints**: Planning constraints
  - `avoid_internal_cuts`: Avoid cutting inside CDS features
  - `dephosphorylate_backbone`: Dephosphorylate vector
  - `min_overhang`: Minimum overhang length (bp)

## Examples

### 1. two_enzyme_cloning.json
Classic restriction cloning with two enzymes for directional insertion.

**Usage:**
```bash
python sim.py --plan-cloning examples/two_enzyme_cloning.json --max-steps 3
```

### 2. golden_gate_assembly.json
Golden Gate (Type IIS) assembly with multiple parts and designed overhangs.

**Usage:**
```bash
python sim.py --plan-cloning examples/golden_gate_assembly.json --prefer-typeIIS --max-steps 2
```

### 3. frame_preserving_fusion.json
Frame-preserving fusion of two fluorescent proteins with a flexible linker.

**Usage:**
```bash
python sim.py --plan-cloning examples/frame_preserving_fusion.json --frame-check --max-steps 4
```

### 4. simple_insertion.yaml
Simple single-gene insertion (YAML format).

**Usage:**
```bash
python sim.py --plan-cloning examples/simple_insertion.yaml --max-steps 2
```

## Advanced Options

### Enzyme Constraints
```bash
# Avoid specific enzymes
python sim.py --plan-cloning spec.json --avoid-enzymes "EcoRI,NotI"

# Restrict to specific enzymes
python sim.py --plan-cloning spec.json --allow-enzymes "BsaI,BsmBI,SpeI,XbaI"
```

### Export Plan
```bash
# Export detailed protocol and intermediate files
python sim.py --plan-cloning spec.json --export-plan output_dir/
```

### Visualizations
```bash
# Include gel simulations for each step
python sim.py --plan-cloning spec.json --print-gels
```

## Creating Your Own Specs

1. Start with one of the examples
2. Update paths to your FASTA files
3. Specify desired assembly order
4. Define junction requirements (directional, frame-preserving, etc.)
5. Set constraints based on your needs

## Tips

- For complex assemblies, increase `--max-steps`
- Use `--prefer-typeIIS` for multi-part assemblies
- Enable `--frame-check` for protein fusions
- Use `--avoid-enzymes` if certain sites are forbidden
- Export plans with `--export-plan` for lab protocols

