# Restriction Enzyme Simulator

A simple Python tool for simulating restriction enzyme cutting on linear DNA sequences. This is Phase 2 of the simulator, now supporting **multiple enzymes** for combined digest analysis, while remaining beginner-friendly with clear explanations and comprehensive comments.

## Features

- **Multiple Enzyme Support**: Now supports 1 to N enzymes for combined digest analysis
- **Extensive Enzyme Database**: Loads from `enzymes.json` file with 160+ enzymes, falls back to 5 built-in enzymes
- **Flexible Input**: Accept DNA sequences as direct strings or from FASTA/text files
- **Case-Insensitive Matching**: Automatically handles both uppercase and lowercase enzyme names and DNA sequences
- **Linear DNA Support**: Calculates fragment lengths for linear DNA molecules
- **Combined Digest Analysis**: Merges cut positions from multiple enzymes and deduplicates overlapping cuts
- **Smart Error Handling**: Helpful suggestions for unknown enzyme names with similarity matching
- **Clear Output**: Detailed analysis showing individual enzyme results and combined digest summary

## Installation

No additional dependencies required! This simulator uses only Python standard library modules:
- `argparse` for command-line arguments
- `json` for loading enzyme database
- `sys` for system operations
- `typing` for type hints

## Enzyme Database

The simulator automatically loads enzymes from `enzymes.json` if present, otherwise uses a built-in database with 5 common enzymes.

### enzymes.json Format

The `enzymes.json` file should contain a JSON array of enzyme objects with the following format:

```json
[
  {
    "name": "EcoRI",
    "site": "GAATTC",
    "cut_index": 1
  },
  {
    "name": "BamHI",
    "site": "GGATCC",
    "cut_index": 1
  }
]
```

**Field descriptions:**
- `name`: The enzyme name (string)
- `site`: The recognition sequence (DNA bases: A, T, C, G)
- `cut_index`: 0-based position where the enzyme cuts within the recognition sequence

**Example:** EcoRI recognizes "GAATTC" and cuts after position 1 (between G and A), so `cut_index` is 1.

## Usage

### Basic Syntax

```bash
python sim.py --seq <DNA_SEQUENCE_OR_FILE> --enz <ENZYME_NAME_1> [ENZYME_NAME_2] ... [ENZYME_NAME_N]
```

### Examples

#### Single enzyme (same as before):
```bash
python sim.py --seq "ATCGAATTCGGGATCCAAA" --enz EcoRI
```

#### Multiple enzymes (NEW!):
```bash
python sim.py --seq "ATCGAATTCGGGATCCAAA" --enz EcoRI BamHI
python sim.py --seq "ATCGAATTCGGGATCCAAA" --enz ecori bamhi hindiii  # Case-insensitive!
```

#### Using a FASTA file with multiple enzymes:
```bash
python sim.py --seq sample_dna.fasta --enz BamHI HindIII
```

### Command-Line Arguments

- `--seq`: DNA sequence as a string or path to a FASTA/text file
- `--enz`: One or more enzyme names (case-insensitive, with helpful suggestions for unknown enzymes)

## Supported Enzymes

With the included `enzymes.json` file, over 160 restriction enzymes are available. Here are examples of the built-in enzymes (available even without `enzymes.json`):

| Enzyme | Recognition Sequence | Cut Position | Example |
|--------|---------------------|--------------|---------|
| EcoRI  | GAATTC              | G^AATTC      | G^AATTC |
| BamHI  | GGATCC              | G^GATCC      | G^GATCC |
| HindIII| AAGCTT              | A^AGCTT      | A^AGCTT |
| PstI   | CTGCAG              | CTGCA^G      | CTGCA^G |
| NotI   | GCGGCCGC            | GC^GGCCGC    | GC^GGCCGC |

*Note: ^ indicates the cut position*

**Additional enzymes available in `enzymes.json`:**
- AluI, ClaI, DraI, EcoRV, HaeIII, HpaI, KpnI, MboI, NcoI, PvuI, SacI, SalI, SmaI, SpeI, SphI, TaqI, XbaI, XhoI, and many more!

To see all available enzymes, run the script with an invalid enzyme name - it will display the complete list.

## Input File Format

The simulator accepts FASTA files with the following format:

```
>Sequence description (optional)
ATCGATCGATCGATCG
```

Or simple text files containing only the DNA sequence:

```
ATCGATCGATCGATCG
```

## Output Examples

### Single Enzyme Example
```
Reading DNA sequence from: ATCGAATTCGGGATCCAAA
DNA sequence length: 19 bp
DNA sequence: ATCGAATTCGGGATCCAAA

Enzyme: EcoRI
Recognition sequence: GAATTC
Cut position: after position 1 in recognition sequence
Found cut positions: [4]

COMBINED DIGEST SUMMARY
========================================
Total cuts (unique across enzymes): 1
Cut positions: [4]
Fragment lengths: [4, 15]
Total length verification: 19 bp (expected: 19 bp)
```

### Multiple Enzymes Example
```
Reading DNA sequence from: ATCGAATTCGGGATCCAAA
DNA sequence length: 19 bp
DNA sequence: ATCGAATTCGGGATCCAAA

Enzyme: EcoRI
Recognition sequence: GAATTC
Cut position: after position 1 in recognition sequence
Found cut positions: [4]

Enzyme: BamHI
Recognition sequence: GGATCC
Cut position: after position 1 in recognition sequence
Found cut positions: [11]

COMBINED DIGEST SUMMARY
========================================
Total cuts (unique across enzymes): 2
Cut positions: [4, 11]
Fragment lengths: [4, 7, 8]
Total length verification: 19 bp (expected: 19 bp)
```

## Error Handling

The simulator handles various error conditions:

- **Invalid enzyme**: Shows similar enzyme suggestions and available enzymes if an unsupported enzyme is specified
- **Case-insensitive matching**: Automatically matches enzyme names regardless of case
- **File not found**: Clear error message if the specified file doesn't exist
- **Invalid DNA characters**: Validates that only A, T, C, G characters are present
- **No cut sites**: Gracefully handles cases where no recognition sequences are found
- **Empty enzyme list**: Provides usage guidance if no enzymes are specified

## File Structure

```
RES/
├── sim.py                 # Main simulator script (Phase 2 - Multi-enzyme support)
├── enzymes.json          # Extended enzyme database (160+ enzymes)
├── sample_dna.fasta      # Sample DNA sequence for testing
├── tests/
│   └── test_multi_enzyme.py  # Test cases for multi-enzyme functionality
├── prompt.txt            # Original requirements
├── prompt 2.txt          # Enhancement requirements
├── prompt 3.txt          # Multi-enzyme requirements
└── README.md             # This file
```

## Testing

Test the simulator with the included sample file:

```bash
# Test with single enzyme
python sim.py --seq sample_dna.fasta --enz EcoRI

# Test with multiple enzymes
python sim.py --seq sample_dna.fasta --enz EcoRI BamHI

# Test with case-insensitive enzyme names
python sim.py --seq "ATCGAATTCGGGATCCAAA" --enz ecori bamhi

# Test with a direct sequence and multiple enzymes
python sim.py --seq "ATCGAATTCGGGATCCAAA" --enz EcoRI BamHI HindIII

# Run automated test suite
python tests/test_multi_enzyme.py
```

## Technical Details

### How It Works

1. **Database Loading**: Loads enzyme database from `enzymes.json` if present, otherwise uses built-in database
2. **Sequence Input**: Reads DNA sequence from string or file, converting to uppercase
3. **Enzyme Validation**: Validates and normalizes enzyme names (case-insensitive) with similarity matching for errors
4. **Individual Analysis**: For each enzyme, finds cut positions using recognition sequences
5. **Cut Position Merging**: Combines cut positions from all enzymes, deduplicating overlapping cuts
6. **Fragment Calculation**: Determines fragment lengths based on merged cut positions for linear DNA
7. **Output Generation**: Displays individual enzyme results and combined digest summary

### Code Structure

- `load_enzyme_database()`: Loads enzyme information from JSON file or returns built-in database
- `read_dna_sequence()`: Handles input from string or file with validation
- `find_cut_sites()`: Locates all recognition sequences in the DNA (original function)
- `find_cut_positions_linear()`: **NEW** - Finds cut positions for a single enzyme
- `merge_cut_positions()`: **NEW** - Combines and deduplicates cut positions from multiple enzymes
- `fragments_linear()`: **NEW** - Calculates fragment lengths for combined digest
- `find_closest_enzyme_names()`: **NEW** - Provides helpful suggestions for unknown enzyme names
- `calculate_fragments()`: Computes fragment sizes for linear DNA (legacy function)
- `main()`: Orchestrates the entire process with multi-enzyme command-line interface

## Future Enhancements (Not in Phase 2)

This is Phase 2 of the simulator. Future phases may include:
- Circular DNA support
- Gel electrophoresis simulation
- Additional enzymes beyond the current 160+
- Graphical output
- Restriction site mapping
- Fragment analysis with sticky ends

## License

This project is open source and available for educational use.

## Contributing

This is a learning project designed for beginners. Feel free to use and modify for educational purposes!
