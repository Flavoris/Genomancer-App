# Restriction Enzyme Simulator

A simple Python tool for simulating restriction enzyme cutting on linear DNA sequences. This is Phase 3 of the simulator, now supporting **advanced IUPAC expansion** and **overhang type classification**, while remaining beginner-friendly with clear explanations and comprehensive comments.

## Features

- **Multiple Enzyme Support**: Supports 1 to N enzymes for combined digest analysis
- **Advanced IUPAC Expansion**: Full support for all IUPAC degenerate bases (A,C,G,T,R,Y,W,S,M,K,B,D,H,V,N)
- **Overhang Type Classification**: Displays overhang type (5' overhang, 3' overhang, Blunt, Unknown) sourced from CSV data
- **Overlapping Matches**: Uses lookahead regex to find all overlapping recognition sites
- **Extensive Enzyme Database**: Loads from `enzymes.json` file with 350+ enzymes, falls back to 5 built-in enzymes
- **Flexible Input**: Accept DNA sequences as direct strings or from FASTA/text files
- **Case-Insensitive Matching**: Automatically handles both uppercase and lowercase enzyme names and DNA sequences
- **Linear DNA Support**: Calculates fragment lengths for linear DNA molecules
- **Combined Digest Analysis**: Merges cut positions from multiple enzymes and deduplicates overlapping cuts
- **Smart Error Handling**: Helpful suggestions for unknown enzyme names with similarity matching
- **Clear Output**: Detailed analysis showing individual enzyme results, overhang information, and combined digest summary

## Installation

No additional dependencies required! This simulator uses only Python standard library modules:
- `argparse` for command-line arguments
- `json` for loading enzyme database
- `re` for pattern matching and IUPAC expansion
- `sys` for system operations
- `unicodedata` for name normalization
- `typing` for type hints

## Enzyme Database

The simulator automatically loads enzymes from `enzymes.json` if present, otherwise uses a built-in database with 5 common enzymes.

### enzymes.json Format

The `enzymes.json` file should contain a JSON array of enzyme objects with the following format:

```json
[
  {
    "name": "AatII",
    "site": "GACGTC",
    "cut_index": 5,
    "overhang_type": "3' overhang"
  }
]
```

**Field descriptions:**
- `name`: The enzyme name (string)
- `site`: The recognition sequence (supports IUPAC degenerate bases: A,C,G,T,R,Y,W,S,M,K,B,D,H,V,N)
- `cut_index`: 0-based position where the enzyme cuts within the recognition sequence
- `overhang_type`: The type of overhang produced ("5' overhang", "3' overhang", "Blunt", or "Unknown")

**Note:** Overhang Type is imported directly from Restriction Enzymes.csv and displayed verbatim (normalized) in output.

**Examples:**
- AatII recognizes "GACGTC" and cuts after position 5, producing a 3' overhang
- EcoRI recognizes "GAATTC" and cuts after position 1 (between G and A), producing a 5' overhang
- PstI recognizes "CTGCAG" and cuts after position 5, producing a 3' overhang

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

#### Using IUPAC degenerate bases:
```bash
python sim.py --seq "ATCGCCACGTGGCCATCG" --enz BsrI
```

#### Using enzymes with overhang information:
```bash
python sim.py --seq "ATCGGCCWGGATCG" --enz EcoRII
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
Site:   GAATTC
Cut @:  index 1
Overhang: 5' overhang
Matches at positions: 4

Fragments:
  [0 - 4] = 4 bp
  [4 - 19] = 15 bp
```

### Multiple Enzymes Example
```
Reading DNA sequence from: ATCGAATTCGGGATCCAAA
DNA sequence length: 19 bp
DNA sequence: ATCGAATTCGGGATCCAAA

Enzyme: EcoRI
Site:   GAATTC
Cut @:  index 1
Overhang: 5' overhang
Matches at positions: 4

Enzyme: BamHI
Site:   GGATCC
Cut @:  index 1
Overhang: 5' overhang
Matches at positions: 11

Fragments:
  [0 - 4] = 4 bp
  [4 - 11] = 7 bp
  [11 - 19] = 8 bp
```

### IUPAC Degenerate Bases Example
```
Reading DNA sequence from: ATCGCCACGTGGCCATCG
DNA sequence length: 18 bp
DNA sequence: ATCGCCACGTGGCCATCG

Enzyme: BsrI
Site:   ACTGGN
Cut @:  index 5
Overhang: Blunt
Matches at positions: 14

Fragments:
  [0 - 14] = 14 bp
  [14 - 18] = 4 bp
```

### Overhang Type Examples
```
Reading DNA sequence from: AAGACGTCTT
DNA sequence length: 10 bp
DNA sequence: AAGACGTCTT

Enzyme: AatII
Site:   GACGTC
Cut @:  index 5
Overhang: 3' overhang
Matches at positions: 2

Fragments:
  [0 - 2] = 2 bp
  [2 - 10] = 8 bp
```

## Error Handling

The simulator handles various error conditions:

- **Invalid enzyme**: Shows similar enzyme suggestions and available enzymes if an unsupported enzyme is specified
- **Case-insensitive matching**: Automatically matches enzyme names regardless of case
- **Ambiguous enzyme names**: Shows all variants when multiple enzymes have the same base name (e.g., EcoRI, EcoRI#2)
- **Duplicate enzyme handling**: Automatically suffixes duplicate enzyme names with #2, #3, etc.
- **File not found**: Clear error message if the specified file doesn't exist
- **Invalid DNA characters**: Validates that only A, T, C, G characters are present
- **Invalid IUPAC characters**: Validates that only supported IUPAC degenerate bases are used
- **No cut sites**: Gracefully handles cases where no recognition sequences are found
- **Empty enzyme list**: Provides usage guidance if no enzymes are specified

## File Structure

```
RES/
├── sim.py                    # Main simulator script (Phase 3 - Enhanced features)
├── enzymes.json              # Extended enzyme database (350+ enzymes)
├── sample_dna.fasta          # Sample DNA sequence for testing
├── tests/
│   ├── test_multi_enzyme.py      # Test cases for multi-enzyme functionality
│   ├── test_enhanced_features.py # Test cases for IUPAC and enhanced features
│   └── test_iupac.py            # Legacy IUPAC tests
├── prompts/
│   ├── prompt.txt               # Original requirements
│   ├── prompt 2.txt             # Enhancement requirements
│   ├── prompt 3.txt             # Multi-enzyme requirements
│   ├── prompt 4.txt             # Advanced features requirements
│   └── prompt 5                 # Current implementation requirements
└── README.md                    # This file
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
```

## Technical Details

### How It Works

1. **Database Loading**: Loads enzyme database from `enzymes.json` if present, otherwise uses built-in database
2. **Sequence Input**: Reads DNA sequence from string or file, converting to uppercase
3. **Enzyme Validation**: Validates and normalizes enzyme names (case-insensitive) with similarity matching for errors
4. **IUPAC Expansion**: Converts IUPAC degenerate bases to regex character classes for pattern matching
5. **Pattern Matching**: Uses lookahead regex to find all overlapping recognition sites
6. **Cut Position Calculation**: Computes break positions using cut_index specifications
7. **Individual Analysis**: For each enzyme, finds cut positions and calculates overhang information
8. **Cut Position Merging**: Combines cut positions from all enzymes, deduplicating overlapping cuts
9. **Fragment Calculation**: Determines fragment lengths based on merged cut positions for linear DNA
10. **Output Generation**: Displays individual enzyme results with overhang information and combined digest summary

### Code Structure

- `iupac_to_regex()`: Converts IUPAC degenerate bases to regex character classes
- `normalize()`: Normalizes enzyme names by removing diacritics and converting to lowercase
- `load_enzyme_database()`: Loads enzyme information from JSON file with support for both legacy and new formats
- `read_dna_sequence()`: Handles input from string or file with validation
- `find_cut_sites()`: Locates all recognition sequences using IUPAC expansion and returns break positions with overhang info
- `find_cut_positions_linear()`: Finds cut positions for a single enzyme
- `merge_cut_positions()`: Combines and deduplicates cut positions from multiple enzymes
- `fragments_linear()`: Calculates fragment lengths for combined digest
- `find_closest_enzyme_names()`: Provides helpful suggestions for unknown enzyme names
- `calculate_fragments()`: Computes fragment sizes for linear DNA (legacy function)
- `main()`: Orchestrates the entire process with enhanced command-line interface

## Future Enhancements (Not in Phase 3)

This is Phase 3 of the simulator. Future phases may include:
- Circular DNA support
- Gel electrophoresis simulation
- Graphical output
- Restriction site mapping
- Fragment analysis with sticky ends
- Support for more complex enzyme behaviors (e.g., methylation sensitivity)

## License

This project is open source and available for educational use.

## Contributing

This is a learning project designed for beginners. Feel free to use and modify for educational purposes!
