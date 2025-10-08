# Restriction Enzyme Simulator

A simple Python tool for simulating restriction enzyme cutting on linear DNA sequences. This is Phase 1 of the simulator, designed to be beginner-friendly with clear explanations and comprehensive comments.

## Features

- **Extensive Enzyme Database**: Loads from `enzymes.json` file with 100+ enzymes, falls back to 5 built-in enzymes
- **Flexible Input**: Accept DNA sequences as direct strings or from FASTA/text files
- **Case-Insensitive**: Automatically handles both uppercase and lowercase DNA sequences
- **Linear DNA Support**: Calculates fragment lengths for linear DNA molecules
- **Clear Output**: Detailed analysis showing cut sites and fragment information
- **Error Handling**: Graceful handling of invalid inputs and missing files

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
python sim.py --seq <DNA_SEQUENCE_OR_FILE> --enz <ENZYME_NAME>
```

### Examples

#### Using a direct DNA sequence:
```bash
python sim.py --seq "ATCGAATTCGGGATCCAAA" --enz EcoRI
```

#### Using a FASTA file:
```bash
python sim.py --seq sample_dna.fasta --enz BamHI
```

### Command-Line Arguments

- `--seq`: DNA sequence as a string or path to a FASTA/text file
- `--enz`: Enzyme name (case-sensitive, available enzymes shown in error message if not found)

## Supported Enzymes

With the included `enzymes.json` file, over 100 restriction enzymes are available. Here are examples of the built-in enzymes (available even without `enzymes.json`):

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

## Output Example

```
Reading DNA sequence from: ATCGAATTCGGGATCCAAA
DNA sequence length: 19 bp
DNA sequence: ATCGAATTCGGGATCCAAA

Enzyme: EcoRI
Recognition sequence: GAATTC
Cut position: after position 1 in recognition sequence

Found 1 cut site(s):
  Cut site 1: position 4

Fragment analysis:
  Number of fragments: 2
  Fragment 1: 5 bp
  Fragment 2: 14 bp

Total length: 19 bp
```

## Error Handling

The simulator handles various error conditions:

- **Invalid enzyme**: Shows available enzymes if an unsupported enzyme is specified
- **File not found**: Clear error message if the specified file doesn't exist
- **Invalid DNA characters**: Validates that only A, T, C, G characters are present
- **No cut sites**: Gracefully handles cases where no recognition sequences are found

## File Structure

```
RES/
├── sim.py                 # Main simulator script
├── enzymes.json          # Extended enzyme database (100+ enzymes)
├── sample_dna.fasta      # Sample DNA sequence for testing
├── prompt.txt            # Original requirements
├── prompt 2.txt          # Enhancement requirements
└── README.md             # This file
```

## Testing

Test the simulator with the included sample file:

```bash
# Test with EcoRI
python sim.py --seq sample_dna.fasta --enz EcoRI

# Test with BamHI
python sim.py --seq sample_dna.fasta --enz BamHI

# Test with a direct sequence
python sim.py --seq "ATCGATCGATCGATCG" --enz HindIII
```

## Technical Details

### How It Works

1. **Database Loading**: Loads enzyme database from `enzymes.json` if present, otherwise uses built-in database
2. **Sequence Input**: Reads DNA sequence from string or file, converting to uppercase
3. **Enzyme Lookup**: Retrieves recognition sequence and cut position from loaded database
4. **Cut Site Detection**: Searches for all occurrences of the recognition sequence (case-insensitive)
5. **Fragment Calculation**: Determines fragment lengths based on cut positions for linear DNA
6. **Output Generation**: Displays detailed analysis including cut sites and fragment information

### Code Structure

- `load_enzyme_database()`: Loads enzyme information from JSON file or returns built-in database
- `read_dna_sequence()`: Handles input from string or file with validation
- `find_cut_sites()`: Locates all recognition sequences in the DNA
- `calculate_fragments()`: Computes fragment sizes for linear DNA
- `main()`: Orchestrates the entire process with command-line interface

## Future Enhancements (Not in Phase 1)

This is Phase 1 of the simulator. Future phases may include:
- Circular DNA support
- Multiple enzyme digestion
- Gel electrophoresis simulation
- Additional enzymes
- Graphical output

## License

This project is open source and available for educational use.

## Contributing

This is a learning project designed for beginners. Feel free to use and modify for educational purposes!
