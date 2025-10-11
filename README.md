# Restriction Enzyme Simulator

A comprehensive Python tool for simulating restriction enzyme cutting on **linear and circular** DNA sequences. Now supporting **circular plasmid digestion**, **advanced IUPAC expansion**, and **overhang type classification**, while remaining beginner-friendly with clear explanations and comprehensive comments.

## Features

- **Linear and Circular DNA Support**: Full support for both linear and circular DNA topology with wrap-around fragment calculation
- **Fragment Sequence Extraction**: Returns actual DNA sequences for each fragment with detailed overhang analysis and FASTA export
- **GenBank & CSV Export (NEW!)**: Export digest results to GenBank format or CSV tables for use in other tools
- **Ligation Compatibility Analysis (NEW!)**: Analyzes which fragment ends are compatible for ligation, with directionality and heuristics
- **Multiple Enzyme Support**: Supports 1 to N enzymes for combined digest analysis
- **Advanced IUPAC Expansion**: Full support for all IUPAC degenerate bases (A,C,G,T,R,Y,W,S,M,K,B,D,H,V,N)
- **Overhang Type Classification**: Displays overhang type (5' overhang, 3' overhang, Blunt, Unknown) from enzyme metadata
- **End Base Analysis**: Calculates and reports the exact bases present at each fragment end based on overhang type
- **Boundary Annotations**: Tracks which enzymes create each fragment boundary with detailed metadata
- **Overlapping Matches**: Uses lookahead regex to find all overlapping recognition sites
- **Extensive Enzyme Database**: Loads from `enzymes.json` file with 350+ enzymes, falls back to 5 built-in enzymes
- **Flexible Input**: Accept DNA sequences as direct strings or from FASTA/text files
- **Case-Insensitive Matching**: Automatically handles both uppercase and lowercase enzyme names and DNA sequences
- **Combined Digest Analysis**: Merges cut positions from multiple enzymes and deduplicates overlapping cuts
- **Smart Error Handling**: Helpful suggestions for unknown enzyme names with similarity matching
- **Clear Output**: Detailed analysis showing fragment positions, wrapping status, boundary enzymes, and validation
- **Restriction Map Visualization**: Text-mode restriction maps showing cut sites along the sequence
- **Gel Simulation**: ASCII agarose gel electrophoresis with multi-lane support, ladders, and circular DNA topology rendering
- **Graphics Output**: Publication-ready SVG/PNG generation for plasmid maps, linear maps, and fragment diagrams

## Installation

**No additional dependencies required!** This simulator uses only Python standard library modules:
- `argparse` for command-line arguments
- `json` for loading enzyme database
- `re` for pattern matching and IUPAC expansion
- `sys` for system operations
- `unicodedata` for name normalization
- `typing` for type hints

**Optional:** For PNG export from SVG graphics, you can install `cairo` and `cairosvg` (see Graphics Output section below).

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
# Linear DNA (default)
python sim.py --seq <DNA_SEQUENCE_OR_FILE> --enz <ENZYME_NAME_1> [ENZYME_NAME_2] ...

# Circular DNA
python sim.py --seq <DNA_SEQUENCE_OR_FILE> --enz <ENZYME_NAME_1> [ENZYME_NAME_2] ... --circular

# Circular DNA with single-cut linearization
python sim.py --seq <DNA_SEQUENCE_OR_FILE> --enz <ENZYME_NAME_1> --circular --circular_single_cut_linearizes
```

### Examples

#### Linear DNA Digestion

##### Single enzyme:
```bash
python sim.py --seq "ATCGAATTCGGGATCCAAA" --enz EcoRI
```

##### Multiple enzymes:
```bash
python sim.py --seq "ATCGAATTCGGGATCCAAA" --enz EcoRI BamHI
python sim.py --seq "ATCGAATTCGGGATCCAAA" --enz ecori bamhi hindiii  # Case-insensitive!
```

##### Using a FASTA file with multiple enzymes:
```bash
python sim.py --seq sample_dna.fasta --enz BamHI HindIII
```

#### Circular DNA Digestion (NEW!)

##### Circular plasmid with two cuts:
```bash
python sim.py --seq plasmid.fasta --enz EcoRI BamHI --circular
```

##### Circular plasmid with single cut (default behavior - intact circle):
```bash
python sim.py --seq plasmid.fasta --enz EcoRI --circular
```

##### Circular plasmid with single cut (linearization):
```bash
python sim.py --seq plasmid.fasta --enz EcoRI --circular --circular_single_cut_linearizes
```

#### Using IUPAC degenerate bases:
```bash
python sim.py --seq "ATCGCCACGTGGCCATCG" --enz BsrI
```

### Command-Line Arguments

#### Basic Arguments

- `--seq`: DNA sequence as a string or path to a FASTA/text file
- `--enz`: One or more enzyme names (case-insensitive, with helpful suggestions for unknown enzymes)
- `--circular`: Treat DNA as circular (default: linear)
- `--circular_single_cut_linearizes`: In circular mode, one cut yields two fragments instead of one intact circle (default: False)

#### Sequence Extraction Arguments (NEW!)

- `--include-seqs`: Include DNA sequences in text output
- `--seq-context <int>`: Limit sequence display to N bases at each end (0 = show full sequence, default: 0)
- `--fasta-out <path>`: Output fragment sequences to a FASTA file

#### Restriction Map Arguments

- `--print-map`: Print restriction map after digestion results
- `--print-map-only`: Print only the restriction map (skip fragment table)
- `--map-width <int>`: Width of the restriction map in characters (default: 80)
- `--map-ticks <int>`: Number of tick marks on the map scale (default: 10)
- `--map-min-hits <int>`: Minimum number of cuts to show an enzyme (default: 1)
- `--map-group-by {enzyme,position}`: Group cuts by enzyme or by position (default: enzyme)
- `--map-show-overhangs`: Show overhang type labels in the map
- `--map-show-sites`: Show recognition sequences in the map
- `--map-circular-origin <int>`: Origin position for circular DNA map (default: 0)

## Fragment Sequence Extraction (NEW!)

The simulator now extracts and reports the actual DNA sequence for each fragment, including detailed information about overhang types and end bases. This feature supports both linear and circular DNA with wrap-around fragments.

### Features

- **Complete Sequence Extraction**: Returns the exact 5'→3' DNA sequence for each fragment
- **Overhang Analysis**: Calculates and reports overhang types (5', 3', or blunt) and overhang lengths
- **End Base Detection**: Shows the specific bases present at each fragment end
- **Circular Wrap-Around**: Correctly handles fragments that wrap around the origin in circular DNA
- **FASTA Export**: Exports fragment sequences to FASTA format with informative headers
- **Context Control**: Option to show only N bases at each end for long fragments

### Basic Usage

#### Display Sequences in Console

```bash
# Show full sequences inline
python sim.py --seq plasmid.fasta --enz EcoRI BamHI --include-seqs

# Show only first/last 10 bases of each fragment (elided display)
python sim.py --seq plasmid.fasta --enz EcoRI BamHI --include-seqs --seq-context 10
```

#### Export to FASTA File

```bash
# Export all fragment sequences to a FASTA file
python sim.py --seq plasmid.fasta --enz EcoRI BamHI --fasta-out fragments.fasta

# Combine with other outputs
python sim.py --seq plasmid.fasta --enz EcoRI BamHI --include-seqs --fasta-out fragments.fasta --print-map
```

### Example Output

#### Console Output with Sequences

```bash
$ python sim.py --seq sample.fasta --enz EcoRI --include-seqs --seq-context 10
```

```
================================================================================
FRAGMENT SEQUENCES
================================================================================

# Fragment 1 (length: 1421 bp)
start=1053  end=2474  mode=linear
ends: left=START, right=EcoRI (5' overhang, 4 bp: AATT)
seq: ATCGATCGAT...GATCGATCGA

# Fragment 2 (length: 2356 bp)
start=2474  end=4830  mode=linear
ends: left=EcoRI (5' overhang, 4 bp: AATT), right=END
seq: GAATTCGCTA...GCTAGCTAGC
```

#### FASTA Output Format

The FASTA file contains one entry per fragment with a structured header:

```
>frag_001|len=1421|start=1053|end=2474|left=START|right=EcoRI:5p:4:AATT
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
...

>frag_002|len=2356|start=2474|end=4830|left=EcoRI:5p:4:AATT|right=END
GAATTCGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
...
```

**Header Format:**
- `frag_NNN`: Fragment ID (001, 002, etc.)
- `len=N`: Fragment length in base pairs
- `start=N`: Start position (0-based, inclusive)
- `end=N`: End position (0-based, exclusive)
- `left=INFO`: Left end information (enzyme:overhang:length:bases or START)
- `right=INFO`: Right end information (enzyme:overhang:length:bases or END)
- `wraps=True`: Present only if fragment wraps around origin (circular DNA)

**End Information Format:**
- `Enzyme:5p:4:AATT` = 5' overhang, 4 bp, bases AATT
- `Enzyme:3p:4:CTGC` = 3' overhang, 4 bp, bases CTGC
- `Enzyme:blunt:0` = Blunt end, no overhang
- `START` = Beginning of linear sequence
- `END` = End of linear sequence

### Circular DNA with Wrap-Around

For circular DNA with multiple cuts, wrap-around fragments are correctly handled:

```bash
$ python sim.py --seq plasmid.fasta --enz EcoRI PstI --circular --include-seqs
```

```
# Fragment 3 (length: 891 bp)
start=4530  end=421  mode=circular
ends: left=PstI (3' overhang, 4 bp: CTGC), right=EcoRI (5' overhang, 4 bp: AATT)
seq: CTGCAGATCG...ATCGGAATTC

Note: This fragment wraps around the origin (position 0)
```

### Overhang Type Examples

The simulator calculates overhang lengths and displays end bases for each fragment:

#### 5' Overhang (EcoRI)
```
Recognition site: GAATTC
Cut position: G^AATTC (after position 1)
Overhang: 4 bp, 5' type
Left fragment end: ...XXXG (trailing 4 bp form 5' overhang)
Right fragment end: AATTCXXX... (recessed end)
```

#### 3' Overhang (PstI)
```
Recognition site: CTGCAG
Cut position: CTGCA^G (after position 5)
Overhang: 4 bp, 3' type
Left fragment end: ...CTGCA (recessed end)
Right fragment end: GXXXX... (3' overhang on complementary strand)
```

#### Blunt End (EcoRV)
```
Recognition site: GATATC
Cut position: GAT^ATC (after position 3)
Overhang: 0 bp, blunt
No overhang bases
```

### Type IIS Enzymes

The simulator supports Type IIS enzymes that cut outside their recognition sequence:

```bash
$ python sim.py --seq test.fasta --enz BsaI --include-seqs
```

BsaI recognizes `GGTCTC` and cuts 7 bp downstream, creating fragments with the cut position offset from the recognition site.

### Validation

The simulator verifies that fragment sequences are correct:
- ✓ Fragment lengths sum to total sequence length
- ✓ Sequences concatenate to reconstruct the original DNA (with appropriate overhang handling)
- ✓ Wrap-around fragments correctly span the origin in circular mode

## Restriction Map Visualization

The simulator now includes a powerful text-mode restriction map that visually summarizes cut sites along the DNA sequence. This feature works for both linear and circular DNA and provides an intuitive overview of where enzymes cut.

### Basic Usage

```bash
# Show restriction map for linear DNA
python sim.py --seq sample_dna.fasta --enz EcoRI BamHI --print-map

# Show only the map (skip fragment table)
python sim.py --seq sample_dna.fasta --enz EcoRI BamHI --print-map-only

# Show map for circular DNA with custom width
python sim.py --seq plasmid.fasta --enz EcoRI BamHI --circular --print-map --map-width 60
```

### Map Features

- **Scale line**: Shows base pair positions along the sequence
- **Ruler line**: Visual guide with tick marks
- **All cuts track**: Shows all cut positions, with `*` marking positions where multiple enzymes cut
- **Enzyme tracks**: Individual tracks for each enzyme showing where it cuts (marked with `^`)
- **Annotations**: Optional overhang types and recognition sequences
- **Wrap-around note**: For circular DNA, shows the wrap-around fragment details

### Example Output - Linear DNA

```bash
$ python sim.py --seq sample_dna.fasta --enz EcoRI HindIII --print-map
```

```
================================================================================
RESTRICTION MAP
================================================================================
Length: 4827 bp   Mode: linear
0        482       965       1447      1930      2412      2895      3377      3860      4827
|---------|---------|---------|---------|---------|---------|---------|---------|---------|
All (5):    ^   ^     ^         *               ^             ^
EcoRI (3):  ^         ^                         ^
HindIII(2):     ^               ^
```

### Example Output - Circular DNA

```bash
$ python sim.py --seq plasmid.fasta --enz BsaI BsmBI --circular --print-map --map-show-overhangs
```

```
================================================================================
RESTRICTION MAP
================================================================================
Length: 3000 bp   Mode: circular
0        300       600       900       1200      1500      1800      2100      2400      3000
|---------|---------|---------|---------|---------|---------|---------|---------|---------|
All (4):  ^     ^                        *                     ^
BsaI (2): ^                 ^            [5' overhang]
BsmBI(2):           ^                     ^  [5' overhang]

wrap spans: last_cut→first_cut = 2890→110 (len = 220)
```

### Advanced Options

#### Show Overhang Types

```bash
python sim.py --seq sample.fasta --enz EcoRI PstI --print-map --map-show-overhangs
```

This adds overhang type labels (e.g., `[5' overhang]`, `[3' overhang]`, `[Blunt]`) to each enzyme track.

#### Show Recognition Sites

```bash
python sim.py --seq sample.fasta --enz EcoRI BamHI --print-map --map-show-sites
```

This displays the recognition sequence next to each enzyme name (truncated if longer than 16 bases).

#### Group by Position

```bash
python sim.py --seq sample.fasta --enz EcoRI BamHI --print-map --map-group-by position
```

Instead of showing one track per enzyme, this groups cuts by position and lists all enzymes that cut at each position:

```
All (5):  ^   ^     ^         *               ^
pos 100: EcoRI [5' overhang]
pos 250: BamHI [5' overhang]
pos 500: EcoRI, HindIII [5' overhang]
...
```

#### Filter Sparse Enzymes

```bash
python sim.py --seq sample.fasta --enz EcoRI BamHI HindIII --print-map --map-min-hits 2
```

Only shows enzymes that cut at least 2 times, useful for focusing on frequent cutters.

#### Adjust Map Width

```bash
# Narrow map for terminal windows
python sim.py --seq sample.fasta --enz EcoRI --print-map --map-width 60

# Wide map for detailed view
python sim.py --seq sample.fasta --enz EcoRI --print-map --map-width 120
```

#### Circular DNA with Custom Origin

```bash
python sim.py --seq plasmid.fasta --enz EcoRI BamHI --circular --print-map --map-circular-origin 1000
```

This adjusts the coordinate system so that position 1000 appears at the left edge of the map.

## Gel Simulation (NEW!)

The simulator now includes an ASCII agarose gel electrophoresis visualization that shows how DNA fragments migrate through a gel. This feature supports multiple lanes, different agarose concentrations, circular DNA topology rendering, and customizable gel parameters.

### Features

- **Visual Migration**: Realistic band migration based on fragment size and agarose percentage
- **Multi-Lane Support**: Compare multiple digests side-by-side with different enzyme combinations
- **Ladder Support**: Three built-in ladder presets (100bp, 1kb, broad)
- **Circular DNA Handling**: Special rendering for supercoiled (SC) and open circular (OC) forms
- **Band Merging**: Automatically merges fragments of similar size with intensity indicators
- **Customizable Parameters**: Adjust gel size, agarose concentration, smear effects, and more

### Basic Usage

```bash
# Simple gel with single digest
python sim.py --seq sample_dna.fasta --enz EcoRI HindIII --simulate-gel

# Show only the gel (skip digestion table)
python sim.py --seq sample_dna.fasta --enz EcoRI HindIII --gel-only

# Adjust agarose percentage
python sim.py --seq sample_dna.fasta --enz EcoRI --simulate-gel --gel-percent 1.2
```

### Gel Parameters

#### Basic Options
- `--simulate-gel`: Add gel visualization after digestion results
- `--gel-only`: Show only the gel (skip digestion table and map)
- `--gel-percent <float>`: Agarose concentration, 0.7-3.0% (default: 1.0)
- `--gel-length <int>`: Gel height in rows (default: 24)
- `--gel-width <int>`: Gel width in characters (default: 80)

#### Ladder Options
- `--gel-ladder <name>`: Ladder preset - "100bp", "1kb", "broad" (default: 1kb)

Available ladders:
- **100bp**: 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000 bp
- **1kb**: 250, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000 bp
- **broad**: 100, 200, 500, 1000, 1500, 2000, 3000, 4000, 5000, 7000, 10000, 12000 bp

#### Advanced Options
- `--gel-lane-gap <int>`: Spacing between lanes (default: 3)
- `--gel-merge-threshold <bp>`: Merge bands closer than this size (default: 20)
- `--gel-smear {none,light,heavy}`: Add gel artifacts (default: none)
- `--gel-dye-front <float>`: Dye front position, 0-1 (default: 0.85)
- `--gel-topology {auto,linearized,native}`: How to render circular DNA (default: auto)

#### Multi-Lane Configuration
- `--lanes-config <json>`: JSON file or string defining multiple lanes

### Example Output - Linear DNA

```bash
$ python sim.py --seq sample_dna.fasta --enz EcoRI HindIII --simulate-gel --gel-percent 1.0 --gel-ladder 1kb
```

```
================================================================================
AGAROSE GEL SIMULATION
================================================================================

┏━┓   ┏━┓
              
    ·     
         ·
              
    •     
         •
              
              
              
         ·
              
              
              
              
    ·     
              
              
              
              
~~~~~~~~~~~~~~

Legend:
  Agarose: 1.0%
  Dye front: 85%

  Ladder: 250, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000 bp
  EcoRI+HindIII: 4872, 1512, 324 bp
```

### Example Output - Circular DNA with No Cuts

For circular DNA with zero cuts, the gel shows supercoiled (SC) and open circular (OC) forms:

```bash
$ python sim.py --seq plasmid.fasta --enz NotI --circular --simulate-gel --gel-topology native
```

```
================================================================================
AGAROSE GEL SIMULATION
================================================================================

┏━┓   ┏━┓
              
    ·     
              
         ·   (OC - slower migration)
              
              
         ·   (SC - faster migration)
              
              
~~~~~~~~~~~~~~

Legend:
  Agarose: 1.0%
  Dye front: 85%

  Ladder: 250, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000 bp
  NotI: 3000 (OC), 3000 (SC) bp
```

### Multi-Lane Gels

Create a JSON configuration file or pass JSON directly to compare multiple digests:

**lanes_config.json:**
```json
[
  {
    "label": "Uncut",
    "enzymes": [],
    "circular": true,
    "notes": "Control plasmid"
  },
  {
    "label": "EcoRI",
    "enzymes": ["EcoRI"],
    "circular": true
  },
  {
    "label": "EcoRI+HindIII",
    "enzymes": ["EcoRI", "HindIII"],
    "circular": true
  }
]
```

```bash
$ python sim.py --seq plasmid.fasta --enz EcoRI --lanes-config lanes_config.json --gel-only
```

This produces a gel with ladder plus three sample lanes showing different digest conditions.

### Gel Topology Modes

For circular DNA, three topology modes control how fragments are rendered:

#### Auto (default)
- **0 cuts**: Renders SC/OC forms (native plasmid topology)
- **1 cut**: Renders as linearized band (single band at full plasmid length)
- **2+ cuts**: Renders as normal fragments

#### Linearized
- Always renders fragments as linear DNA, even for single cuts

#### Native
- Preserves circular topology rendering
- Shows SC/OC forms for 0-1 cuts
- Use with `--gel-topology native` flag

Example:
```bash
# Force native topology rendering
python sim.py --seq plasmid.fasta --enz EcoRI --circular --simulate-gel --gel-topology native
```

### Band Intensity and Merging

Fragments that migrate to the same position or are within the merge threshold are combined into a single brighter band:

- `·` - Single fragment (light)
- `•` - Two fragments merged
- `▮` - Three fragments merged
- `█` - Four or more fragments merged

The legend shows the multiplicity, e.g., "500 (3×)" indicates three fragments near 500 bp merged into one band.

### Agarose Percentage Effects

Different agarose concentrations provide optimal resolution for different size ranges:

- **0.7-0.8%**: Best for large fragments (1-20 kb)
- **1.0-1.2%**: General purpose (0.5-10 kb)
- **1.5-2.0%**: Small fragments (0.1-3 kb)
- **2.5-3.0%**: Very small fragments (50-1000 bp)

Example:
```bash
# High percentage gel for small fragments
python sim.py --seq sample.fasta --enz MseI --simulate-gel --gel-percent 2.0 --gel-ladder 100bp
```

### Smear Effects

Simulate gel artifacts and DNA degradation:

```bash
# Light smear (trailing below bands)
python sim.py --seq sample.fasta --enz EcoRI --simulate-gel --gel-smear light

# Heavy smear (more artifacts, overloading)
python sim.py --seq sample.fasta --enz EcoRI --simulate-gel --gel-smear heavy
```

## Export (GenBank & CSV) (NEW!)

The simulator now supports exporting restriction digest results to **GenBank** (.gb/.gbk) and **CSV** formats. This allows you to save digest results for use in other tools or for archival purposes.

### Features

- **GenBank Export**: Standard-compliant GenBank files with LOCUS, FEATURES, and ORIGIN sections
- **CSV Export**: Separate CSV files for fragments and cut sites with full metadata
- **Circular Topology Support**: Handles circular DNA with proper `join()` notation for wrap-around features
- **Complete Metadata**: Includes enzyme names, recognition sites, overhang types, and positions
- **Customizable**: Add custom definitions, organism names, and topology overrides

### GenBank Export

Export your digest results to a GenBank-formatted file that can be opened in tools like SnapGene, Benchling, or ApE:

```bash
# Basic GenBank export (linear DNA)
python sim.py --seq sample.fasta --enz EcoRI BamHI --export-genbank output.gbk

# Circular DNA export
python sim.py --seq plasmid.fasta --enz EcoRI BamHI --circular --export-genbank plasmid_digest.gb

# With custom definition and organism
python sim.py --seq plasmid.fasta --enz EcoRI --circular \
  --export-genbank pUC19_EcoRI.gbk \
  --gb-definition "pUC19 digested with EcoRI" \
  --source "E. coli cloning vector"
```

#### GenBank File Structure

The exported GenBank file includes:

1. **LOCUS**: Sanitized name (max 16 chars), length, DNA type, topology, and date
2. **DEFINITION**: Custom or default description
3. **ACCESSION/VERSION**: Blank (standard for synthetic constructs)
4. **SOURCE**: Organism name with "synthetic construct" lineage
5. **FEATURES**:
   - `source` feature spanning the entire molecule
   - `misc_feature` for each restriction site with enzyme name and overhang details
   - `misc_feature` for each fragment with length and boundary information
6. **ORIGIN**: Formatted sequence (60 nt/line, 10-nt blocks, 1-based indexing)

#### Circular DNA Handling

For circular DNA with wrap-around features, the GenBank export uses `join()` notation:

```
misc_feature    join(901..1000,1..100)
                /label="fragment_1"
                /note="length=200bp; left=EcoRI(5' overhang), right=BamHI(5' overhang)"
```

### CSV Export

Export fragment and cut site data to CSV files for analysis in spreadsheet tools:

```bash
# Basic CSV export (creates two files: prefix_fragments.csv and prefix_cuts.csv)
python sim.py --seq sample.fasta --enz EcoRI BamHI --export-csv results/digest

# This creates:
#   results/digest_fragments.csv
#   results/digest_cuts.csv

# Circular DNA export
python sim.py --seq plasmid.fasta --enz EcoRI --circular --export-csv plasmid_analysis
```

#### CSV File Formats

**Fragments CSV** (`<prefix>_fragments.csv`):
- `fragment_id`: Fragment number
- `start_idx`: Start position (0-based, inclusive)
- `end_idx`: End position (0-based, exclusive)
- `mode`: "linear" or "circular"
- `length`: Fragment length in base pairs
- `left_enzyme`: Enzyme at left boundary (empty if none)
- `left_overhang_type`: Type of left overhang
- `left_overhang_len`: Length of left overhang
- `left_end_bases`: Sequence of left overhang bases
- `right_enzyme`: Enzyme at right boundary (empty if none)
- `right_overhang_type`: Type of right overhang
- `right_overhang_len`: Length of right overhang
- `right_end_bases`: Sequence of right overhang bases
- `sequence`: Full fragment sequence

**Cuts CSV** (`<prefix>_cuts.csv`):
- `cut_id`: Sequential cut number
- `pos`: Cut position (0-based)
- `enzyme`: Enzyme name
- `recognition_site`: Recognition sequence
- `cut_index`: Position within site where cut occurs
- `overhang_type`: "5' overhang", "3' overhang", or "Blunt"
- `overhang_len`: Length of overhang in base pairs

### Export Options

All export-related command-line flags:

```bash
--export-genbank <path>            # Export to GenBank file at specified path
--export-csv <prefix>              # Export to CSV files (creates <prefix>_fragments.csv and <prefix>_cuts.csv)
--gb-definition "<text>"           # GenBank DEFINITION line (default: "Restriction digest export")
--source "<organism>"              # GenBank SOURCE organism (default: "synthetic DNA")
--topology {linear,circular}       # Override topology for export (defaults to current run mode)
```

### Combined Export Example

You can export to both formats simultaneously:

```bash
# Export to both GenBank and CSV
python sim.py --seq plasmid.fasta --enz EcoRI BamHI HindIII --circular \
  --export-genbank output/plasmid_digest.gbk \
  --export-csv output/plasmid_digest \
  --gb-definition "Triple digest of pUC19" \
  --source "pUC19 cloning vector"

# This creates:
#   output/plasmid_digest.gbk
#   output/plasmid_digest_fragments.csv
#   output/plasmid_digest_cuts.csv
```

### Important Notes

#### Coordinate Systems

- **GenBank**: Uses 1-based, inclusive-inclusive coordinates (e.g., `1..100`)
- **CSV**: Uses 0-based, inclusive-exclusive coordinates (e.g., `start=0, end=100`)
- **Circular wraps**: GenBank uses `join()` notation; CSV shows actual start/end with `mode=circular`

#### Fragment Sequences

In CSV export, the `sequence` column contains the full fragment sequence:
- Linear fragments: Direct subsequence extraction
- Circular wraps: Concatenated sequence from wrap-around (end of sequence + beginning of sequence)

#### Overhang Information

Both formats include complete overhang metadata:
- Overhang type (5', 3', or Blunt)
- Overhang length in base pairs
- End bases (5'→3' orientation) where applicable

## Ligation Compatibility Analysis (NEW!)

The simulator now includes **ligation compatibility analysis** that determines which fragment ends can be ligated together. This feature analyzes sticky-end compatibility, directionality, and provides ligation heuristics (GC%, Tm) for cloning applications.

### Features

- **Sticky-End Compatibility**: Checks if fragment ends can anneal based on overhang sequence complementarity
- **Directionality Analysis**: Identifies directional vs. palindromic overhangs (important for preventing self-ligation)
- **Blunt-End Support**: Optional analysis of blunt-blunt compatibility
- **Ligation Heuristics**: Calculates GC%, Tm (Wallace rule), and provides ligation guidance
- **Multiple Output Formats**: Pairs list, compatibility matrix, or detailed report
- **JSON Export**: Machine-readable output for downstream analysis

### Compatibility Rules

The simulator applies the following rules for ligation compatibility:

1. **Sticky ↔ Sticky only**: Sticky ends can only ligate with other sticky ends (blunt is optional)
2. **Sticky ↔ Blunt is incompatible**: Mixed sticky/blunt ends cannot ligate
3. **Length must match**: Both overhangs must be the same length (e.g., 4 bp with 4 bp)
4. **Type must match**: Both must be 5' overhangs OR both must be 3' overhangs
5. **Complementarity**: The sticky sequences must be reverse complements: `sticky_seq_A == revcomp(sticky_seq_B)`
6. **Directionality**: A pair is directional if the overhang is non-palindromic (prevents self-ligation)

### Basic Usage

```bash
# Basic compatibility analysis
python sim.py --seq plasmid.fasta --enz EcoRI HindIII --compatibility

# Show detailed information with heuristics
python sim.py --seq plasmid.fasta --enz EcoRI HindIII --compatibility --compat-summary detailed

# Matrix view of all end-to-end compatibility
python sim.py --seq plasmid.fasta --enz EcoRI HindIII --compatibility --compat-summary matrix

# Filter to directional pairs only (non-palindromic)
python sim.py --seq plasmid.fasta --enz SpeI XbaI --compatibility --require-directional

# Include blunt-blunt as compatible
python sim.py --seq plasmid.fasta --enz EcoRI EcoRV --compatibility --include-blunt

# Export to JSON for downstream analysis
python sim.py --seq plasmid.fasta --enz EcoRI HindIII --compatibility --json-out results.json
```

### Compatibility Options

All ligation compatibility flags:

```bash
--compatibility                 # Enable compatibility analysis
--compat-summary {pairs,matrix,detailed}  # Output format (default: pairs)
--require-directional          # Filter to directional pairs only
--include-blunt                # Include blunt-blunt as compatible (default: False)
--min-overhang N               # Minimum overhang length for sticky classification (default: 1)
--json-out <path>              # Export results to JSON file
```

### Output Formats

#### Pairs Format (Default)

Lists all compatible end pairs with ligation details:

```
================================================================================
LIGATION COMPATIBILITY ANALYSIS - COMPATIBLE PAIRS
================================================================================

Compatible pair #1 (k=4):
  [frag1:right] EcoRI 5' overhang: AATT
  ↔ [frag2:left] EcoRI 5' overhang: AATT (revcomp match)
  directionality: NO (palindromic), GC%: 0.0, Tm≈8°C

Compatible pair #2 (k=4):
  [frag3:right] SpeI 5' overhang: CTAG
  ↔ [frag7:left] XbaI 5' overhang: CTAG (revcomp match)
  directionality: NO (palindromic), GC%: 50.0, Tm≈12°C

Total compatible pairs: 2
```

#### Matrix Format

Shows N×N compatibility matrix for all fragment ends:

```
================================================================================
LIGATION COMPATIBILITY ANALYSIS - COMPATIBILITY MATRIX
================================================================================

     F0R  F1L  F1R  F2L
    --------------------
F0R   ·   ✓   .   . 
F1L   ✓   ·   .   . 
F1R   .   .   ·   ✓ 
F2L   .   .   ✓   · 

Legend:
  ✓  = compatible (sticky ends)
  •  = compatible (blunt ends)
  .  = incompatible
  ·  = same end

Total compatible pairs: 2
```

#### Detailed Format

Includes all metadata, heuristics, and scores:

```
================================================================================
LIGATION COMPATIBILITY ANALYSIS - DETAILED REPORT
================================================================================

Pair #1:
  End A: Fragment 1 (right), Enzyme: EcoRI
         Type: 5' overhang, Length: 4 bp
         Sequence: 5'-AATT-3'
         GC%: 0.0%, Tm: 8.0°C
  End B: Fragment 2 (left), Enzyme: EcoRI
         Type: 5' overhang, Length: 4 bp
         Sequence: 5'-AATT-3'
         GC%: 0.0%, Tm: 8.0°C
  Compatible: True
  Directional: False
  Note: 5' overhang, 4 bp overhang, non-directional (palindromic)

Total compatible pairs: 1
```

### Example Scenarios

#### Classic Compatible Pair: SpeI and XbaI

```bash
python sim.py --seq test.fasta --enz SpeI XbaI --compatibility
```

SpeI (`A↓CTAGT`) and XbaI (`T↓CTAGA`) both produce `CTAG` overhangs, making them compatible for creating directional cloning sites.

#### Incompatible: EcoRI and MfeI

```bash
python sim.py --seq test.fasta --enz EcoRI MfeI --compatibility
```

Although EcoRI (`G↓AATTC`) and MfeI (`C↓AATTG`) produce similar overhangs, they may have different lengths or complementarity depending on the exact cut pattern.

#### 3' Overhangs: PstI and NsiI

```bash
python sim.py --seq test.fasta --enz PstI NsiI --compatibility
```

PstI (`CTGCA↓G`) produces a 3' overhang. Compatible ends must also have 3' overhangs with matching sequences.

#### Blunt End Ligation

```bash
python sim.py --seq test.fasta --enz EcoRV SmaI --compatibility --include-blunt
```

Both EcoRV and SmaI are blunt cutters. With `--include-blunt`, their ends are compatible for blunt-end ligation (requires T4 DNA ligase).

### Heuristics Explained

#### GC Percentage

Shows the GC content of the overhang sequence:
- Higher GC% = stronger annealing
- Typical range: 0-100%

#### Tm Estimation (Wallace Rule)

Rough melting temperature for short overhangs:
- Formula: `Tm ≈ 2×(A+T) + 4×(G+C)`
- Only accurate for short oligos (<14 nt)
- Higher Tm = more stable annealing
- Useful for choosing ligation temperatures

Example:
- `AATT`: Tm ≈ 8°C (weak, AT-rich)
- `GGCC`: Tm ≈ 16°C (strong, GC-rich)

### Directionality

**Directional pairs** have non-palindromic overhangs that enforce insert orientation and prevent self-ligation:

```bash
# Filter to directional pairs only
python sim.py --seq test.fasta --enz SpeI XbaI --compatibility --require-directional
```

**Non-directional (palindromic)** overhangs can self-ligate and don't enforce orientation:
- EcoRI: `AATT` (palindrome)
- BamHI: `GATC` (palindrome)

### JSON Export

Export compatibility results for programmatic analysis:

```bash
python sim.py --seq test.fasta --enz EcoRI HindIII --compatibility --json-out compat.json
```

JSON format:
```json
[
  {
    "end_a": {
      "fragment_id": 0,
      "polarity": "right",
      "enzyme": "EcoRI",
      "overhang_type": "5' overhang",
      "overhang_len": 4,
      "sticky_seq": "AATT",
      "gc_percent": 0.0,
      "tm": 8.0,
      "position": 100
    },
    "end_b": {
      "fragment_id": 1,
      "polarity": "left",
      "enzyme": "EcoRI",
      "overhang_type": "5' overhang",
      "overhang_len": 4,
      "sticky_seq": "AATT",
      "gc_percent": 0.0,
      "tm": 8.0,
      "position": 200
    },
    "compatible": true,
    "directional": false,
    "note": "5' overhang, 4 bp overhang, non-directional (palindromic)"
  }
]
```

### Circular DNA Support

Compatibility analysis works with both linear and circular DNA:

```bash
# Circular plasmid with multiple cuts
python sim.py --seq plasmid.fasta --enz EcoRI HindIII --circular --compatibility

# Linearized plasmid (single cut)
python sim.py --seq plasmid.fasta --enz EcoRI --circular --circular_single_cut_linearizes --compatibility
```

For circular DNA with 2+ cuts, wrap-around fragment ends are analyzed correctly.

### Edge Cases

The simulator handles several edge cases:

1. **Mixed 5' and 3' overhangs**: Always incompatible (different geometry)
2. **Unequal lengths**: Incompatible even if partially complementary
3. **Same cut site**: Ends from the same cut can be compatible (useful for re-ligation)
4. **Coincident cuts**: Multiple enzymes at the same position are handled correctly
5. **Palindromic sticky ends**: Compatible but flagged as non-directional

### Practical Applications

#### Subcloning Strategy

```bash
# Check if insert and vector ends are compatible
python sim.py --seq insert.fasta --enz EcoRI XbaI --compatibility
python sim.py --seq vector.fasta --enz EcoRI XbaI --compatibility
```

#### Directional Cloning

```bash
# Find directional enzyme pairs
python sim.py --seq gene.fasta --enz SpeI XbaI --compatibility --require-directional
```

#### Blunt-End Cloning

```bash
# Check blunt-end compatibility
python sim.py --seq pcr_product.fasta --enz SmaI EcoRV --compatibility --include-blunt
```

## Theoretical Enzyme Compatibility (No Digest) (NEW!)

The simulator now supports **theoretical compatibility analysis** that predicts sticky-end compatibility between enzymes without requiring a DNA sequence or digest. This feature analyzes enzyme metadata (recognition sites, cut positions, overhang types) to compute compatibility, making it ideal for planning cloning strategies before running experiments.

### Features

- **No Sequence Required**: Analyze enzyme compatibility without loading or digesting DNA
- **IUPAC Template Support**: Handles degenerate bases in recognition sites (useful for Type IIS enzymes)
- **Overhang Derivation**: Automatically computes sticky-end templates from enzyme metadata
- **Palindrome Detection**: Identifies directional vs. non-directional overhangs
- **Pairwise Analysis**: Compare specific enzymes or generate full compatibility matrix
- **Multiple Output Formats**: Pairs list, matrix view, or detailed JSON output
- **Batch Processing**: Analyze all enzymes in database with `--theoretical-all`

### How It Works

The theoretical mode:

1. Loads enzyme metadata from `enzymes.json` (recognition site, cut_index, overhang_type)
2. Derives the sticky-end template (IUPAC 5'→3' sequence) for each enzyme
3. Calculates overhang length (k) based on cut geometry
4. Determines if each overhang is palindromic (self-complementary)
5. Checks pairwise IUPAC-aware compatibility between all enzyme pairs
6. Filters results based on your criteria (directionality, minimum length, etc.)

### Compatibility Rules (Theoretical Mode)

1. **Sticky ↔ Sticky only**: Blunt-blunt optional with `--include-blunt`
2. **Length must match**: Overhang lengths must be identical
3. **Type must match**: Both 5' OR both 3' overhangs
4. **IUPAC compatibility**: Templates must be complementary considering degenerate bases
5. **Geometry match**: 5' vs 3' geometry must be the same

### Basic Usage

```bash
# Compare specific enzymes
python sim.py --theoretical-enzymes "EcoRI,XbaI,MfeI"

# Generate compatibility matrix for all enzymes
python sim.py --theoretical-all --format matrix

# Filter to directional pairs only
python sim.py --theoretical-enzymes "EcoRI,XbaI,SpeI,NheI" --require-directional

# Detailed JSON output
python sim.py --theoretical-enzymes "EcoRI,MfeI" --format detailed --json-out results.json
```

### Theoretical Mode Options

All theoretical compatibility flags:

```bash
--theoretical-enzymes "Enz1,Enz2,..."  # Comma-separated enzyme names
--theoretical-all                       # Analyze all enzymes in database
--format {pairs,matrix,detailed}        # Output format (default: pairs)
--require-directional                   # Only show non-palindromic pairs
--include-blunt                         # Include blunt-blunt as compatible
--min-overhang N                        # Minimum overhang length (default: 1)
--json-out <path>                       # Export to JSON file
```

### Output Formats

#### Pairs Format (Default)

Shows compatible enzyme pairs with templates and palindromicity:

```
Theoretical sticky-end compatibility (no digest)
================================================================================

EcoRI  | 5' overhang k=4 | tpl=AATT | palindromic: YES
MfeI   | 5' overhang k=4 | tpl=AATT | palindromic: YES
  ✔ Compatible (non-directional)

SpeI   | 5' overhang k=4 | tpl=ACTA | palindromic: NO
XbaI   | 5' overhang k=4 | tpl=CTAG | palindromic: NO
  ✔ Compatible (directional)

Total compatible pairs: 2
```

#### Matrix Format

Shows N×N compatibility grid:

```
Theoretical compatibility matrix (no digest)
================================================================================

       EcoRI  MfeI  SpeI  XbaI  PstI
     -----------------------------------
EcoRI    .     ✓     .     .     .   
MfeI     ✓     .     .     .     .   
SpeI     .     .     .     ▶     .   
XbaI     .     .     ▶     .     .   
PstI     .     .     .     .     .   

Legend:
  ✓ = compatible (non-directional/palindromic)
  ▶ = compatible (directional)
  . = incompatible or same enzyme

Total compatible pairs: 2
```

#### Detailed Format (JSON)

Complete metadata for programmatic analysis:

```json
[
  {
    "enzyme_a": "EcoRI",
    "enzyme_b": "MfeI",
    "overhang_type": "5' overhang",
    "k": 4,
    "template_a": "AATT",
    "template_b": "AATT",
    "compatible": true,
    "directional": false,
    "reason": "Compatible",
    "palindromic_a": true,
    "palindromic_b": true
  }
]
```

### Example Use Cases

#### Planning Directional Cloning

```bash
# Find all directional enzyme pairs in database
python sim.py --theoretical-all --require-directional --format matrix
```

Identifies enzyme pairs that enforce insert orientation (e.g., SpeI/XbaI, NheI/XbaI).

#### Checking Isoschizomer Compatibility

```bash
# Are EcoRI and MfeI compatible?
python sim.py --theoretical-enzymes "EcoRI,MfeI"
```

Output shows they produce identical `AATT` overhangs and are compatible.

#### Surveying Database for Compatible Partners

```bash
# Find all 4bp 5' overhang enzymes
python sim.py --theoretical-all --min-overhang 4 --format matrix
```

Generates a matrix of all enzymes with ≥4bp overhangs.

#### Type IIS Enzyme Analysis

```bash
# Analyze Type IIS enzymes (may have N in templates)
python sim.py --theoretical-enzymes "BsaI,Esp3I" --format detailed
```

IUPAC-aware matching handles degenerate bases in recognition sites.

### IUPAC Support

The theoretical mode fully supports IUPAC degenerate bases:

- **Standard bases**: A, C, G, T
- **Degenerate codes**: R (A/G), Y (C/T), S (G/C), W (A/T), K (G/T), M (A/C), B (C/G/T), D (A/G/T), H (A/C/T), V (A/C/G), N (any)

Two templates are compatible if at every position, their IUPAC sets intersect when accounting for reverse complementarity.

### Theoretical vs. Actual Overhangs

**Important Notes:**

1. **Type II enzymes**: Theoretical templates accurately predict actual overhangs
2. **Type IIS enzymes**: Overhangs extend beyond recognition sites, so templates may contain `N` (unknown bases)
3. **Genomic context**: Type IIS actual sequences depend on flanking DNA, which theoretical mode doesn't have
4. **Use case**: Theoretical mode predicts *potential* compatibility; actual digest confirms *real* compatibility

For Type IIS enzymes or when you need exact sequences, use regular digest mode with `--compatibility`.

### Large Database Warning

When using `--theoretical-all` with large enzyme databases:

```bash
python sim.py --theoretical-all --format matrix
```

Output:
```
Analyzing all 350 enzymes in database...
Warning: Large output expected (350 enzymes)
```

The matrix will be large but remains useful for identifying clusters of compatible enzymes.

### Comparison with Digest Mode

| Feature | Theoretical Mode | Digest Mode |
|---------|-----------------|-------------|
| **Requires sequence** | No | Yes |
| **Speed** | Instant | Depends on sequence length |
| **Use case** | Planning, database survey | Actual cloning, verification |
| **Overhang accuracy** | Template (may have N) | Exact sequences |
| **Type IIS support** | Partial (N in templates) | Full (actual flanking bases) |

**Recommendation**: Use theoretical mode for planning and exploration, then confirm with digest mode before experimental work.

## Graphics Output (SVG/PNG) (NEW!)

The simulator now includes **publication-ready SVG graphics generation** for plasmid maps, linear restriction maps, and fragment diagrams. SVG files work in all modern applications and require **no additional dependencies**. Optional PNG export is available if needed.

### Features

- **Circular Plasmid Maps**: Render plasmids as circular maps with radial cut markers and labels
- **Linear Restriction Maps**: Traditional linear maps with cut positions and enzyme labels
- **Fragment Diagrams**: Visual representation of fragment sizes with wrap-around indication
- **Publication Ready**: Clean, professional styling with customizable themes
- **No Dependencies**: SVG generation works out-of-the-box with no external libraries
- **Deterministic Colors**: Enzyme colors are consistent across runs (hash-based coloring)
- **Collision Avoidance**: Smart label placement to prevent overlapping text
- **Overhang Badges**: Shows overhang type (5', 3', B for blunt) directly from enzymes.json
- **PNG Export (Optional)**: High-DPI PNG conversion available with cairo/cairosvg

### Basic Usage

Generate SVG graphics by adding output path flags (no additional dependencies required):

```bash
# Generate plasmid map
python sim.py --seq plasmid.fasta --enz EcoRI HindIII --circular --out-svg map.svg

# Generate linear restriction map
python sim.py --seq dna.fasta --enz EcoRI BamHI --out-svg-linear linear.svg

# Generate fragment diagram
python sim.py --seq plasmid.fasta --enz EcoRI HindIII --circular --out-svg-fragments frags.svg

# Generate all three graphics at once
python sim.py --seq plasmid.fasta --enz EcoRI HindIII --circular \
  --out-svg map.svg \
  --out-svg-linear linear.svg \
  --out-svg-fragments frags.svg
```

**SVG files can be used directly** in web browsers, vector editors (Inkscape, Illustrator), presentations (PowerPoint, Keynote), documents (Word, Google Docs), and publications. They scale to any size without quality loss.

### Graphics Options

All graphics generation flags:

```bash
--out-svg <path>              # Plasmid/DNA map SVG output path
--out-svg-linear <path>       # Linear restriction map SVG output path
--out-svg-fragments <path>    # Fragment diagram SVG output path
--png                         # Also generate PNG files (optional, requires cairo library)
--theme {light,dark}          # Color theme (default: light)
--title <text>                # Custom title for graphics
--show-sites                  # Include recognition sequences on labels
--hide-overhangs              # Hide overhang type badges (shown by default)
--origin <position>           # Circular map origin position (default: 0)
--svg-width <pixels>          # Override width for linear/fragment diagrams (default: 900)
--svg-height <pixels>         # Override height for linear/fragment diagrams
```

### Examples

#### Complete Example with All Options

```bash
python sim.py --seq plasmid.fasta --enz EcoRI HindIII BamHI --circular \
  --out-svg map.svg --out-svg-linear linear.svg --out-svg-fragments frags.svg \
  --theme dark --show-sites --title "pUC19 Digest" --origin 500
```

This command:
- Digests circular plasmid with three enzymes
- Generates all three types of graphics (SVG format)
- Uses dark theme for better presentation visibility
- Shows recognition sequences on labels
- Uses custom title "pUC19 Digest"
- Sets circular origin at position 500

Add `--png` if you also want PNG files (requires Cairo library to be installed).

#### Dark Theme with Custom Dimensions

```bash
python sim.py --seq dna.fasta --enz EcoRI --out-svg-linear linear.svg \
  --theme dark --svg-width 1200 --svg-height 200
```

#### Show Recognition Sites and Overhang Types

```bash
python sim.py --seq plasmid.fasta --enz EcoRI PstI NotI --circular \
  --out-svg map.svg --show-sites
```

The labels will show:
- Enzyme names (e.g., "EcoRI")
- Recognition sequences (e.g., "(GAATTC)") when `--show-sites` is used
- Overhang badges (e.g., "[5']", "[3']", "[B]") by default
- Position (e.g., "@ 1234")

### PNG Export (Optional)

**SVG files are the recommended output format** - they work in all modern applications, scale perfectly, and require no additional dependencies. However, if you specifically need PNG format, you can optionally enable PNG export.

To generate high-resolution PNG images alongside SVG files:

1. Install the Cairo graphics library on your system:
   ```bash
   # macOS
   brew install cairo
   
   # Ubuntu/Debian
   sudo apt-get install libcairo2-dev
   
   # Windows
   # Download from https://www.cairographics.org/
   ```

2. Install the `cairosvg` Python package:
   ```bash
   pip install cairosvg
   ```

3. Use the `--png` flag:
   ```bash
   python sim.py --seq plasmid.fasta --enz EcoRI --circular \
     --out-svg map.svg --png
   ```

This creates both `map.svg` and `map.png` (at 2× resolution for high DPI displays).

**If Cairo is not installed:** The simulator will still save the SVG file successfully and display a helpful warning message. Your graphics are not affected - you can use the SVG files directly in publications, presentations, and web pages.

### Plasmid Map Details

Circular plasmid maps show:
- **Main circle**: Represents the plasmid backbone
- **Tick marks**: Major ticks every 1000 bp (labeled), minor ticks at appropriate intervals
- **Cut markers**: Radial lines from the circle indicating restriction sites
- **Enzyme labels**: Positioned around the circle with collision avoidance
- **Position labels**: Show cut positions (e.g., "EcoRI @ 234")
- **Multiple enzymes at same position**: Indicated with "×N" marker
- **Origin handling**: Use `--origin` to rotate the map (useful for placing features at top)

### Linear Map Details

Linear restriction maps show:
- **Horizontal ruler**: Scale bar with position markers
- **Triangular markers**: Point to cut positions on the ruler (▼)
- **Enzyme labels**: Alternate above/below ruler to avoid collisions
- **Multi-enzyme sites**: Shows all enzymes when multiple cut at same position

### Fragment Diagram Details

Fragment diagrams show:
- **Scaled blocks**: Rectangles sized proportionally to fragment length
- **Size labels**: Fragment sizes in bp (e.g., "1,234 bp")
- **Wrap fragments**: Split blocks with dashed connector for circular DNA wrap-around
- **Color coding**: Different colors for normal vs. wrap fragments
- **Scale ruler**: Position reference at bottom

### Theme Comparison

**Light theme** (default):
- White background
- Black text and borders
- Muted colors for enzyme markers

**Dark theme**:
- Near-black background (#1a1a1a)
- Light text (#e0e0e0)
- High-contrast colors for better visibility

Example:
```bash
# Light theme (default)
python sim.py --seq dna.fasta --enz EcoRI --out-svg map_light.svg

# Dark theme
python sim.py --seq dna.fasta --enz EcoRI --out-svg map_dark.svg --theme dark
```

### Integration with Other Features

Graphics output works seamlessly with other simulator features:

```bash
# Combine with text restriction map and gel simulation
python sim.py --seq plasmid.fasta --enz EcoRI HindIII --circular \
  --print-map --simulate-gel \
  --out-svg map.svg --out-svg-fragments frags.svg --png
```

This command produces:
1. Text-based restriction map (console output)
2. ASCII gel simulation (console output)
3. SVG plasmid map (map.svg)
4. PNG plasmid map (map.png)
5. SVG fragment diagram (frags.svg)
6. PNG fragment diagram (frags.png)

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

### Linear DNA - Single Enzyme Example
```
Reading DNA sequence from: ATCGAATTCGGGATCCAAA
Topology: linear
DNA sequence length: 19 bp
DNA sequence: ATCGAATTCGGGATCCAAA

Enzyme: EcoRI
Site:   GAATTC
Cut @:  index 1
Overhang: 5' overhang
Matches at positions: 4

================================================================================
DIGESTION RESULTS
================================================================================
Mode: linear
Sequence length: 19 bp
Total cuts: 1
Cut positions: 4
Fragments generated: 2

Fragment Details:
--------------------------------------------------------------------------------
Index    Start    End      Length     Wraps    Boundaries
--------------------------------------------------------------------------------
0        0        4        4          No       START -> 4(EcoRI)
1        4        19       15         No       4(EcoRI) -> END
--------------------------------------------------------------------------------
✓ Fragment lengths sum correctly to 19 bp
```

### Linear DNA - Multiple Enzymes Example
```
Reading DNA sequence from: ATCGAATTCGGGATCCAAA
Topology: linear
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

================================================================================
DIGESTION RESULTS
================================================================================
Mode: linear
Sequence length: 19 bp
Total cuts: 2
Cut positions: 4, 11
Fragments generated: 3

Fragment Details:
--------------------------------------------------------------------------------
Index    Start    End      Length     Wraps    Boundaries
--------------------------------------------------------------------------------
0        0        4        4          No       START -> 4(EcoRI)
1        4        11       7          No       4(EcoRI) -> 11(BamHI)
2        11       19       8          No       11(BamHI) -> END
--------------------------------------------------------------------------------
✓ Fragment lengths sum correctly to 19 bp
```

### Circular DNA - Two Cuts with Wrap-Around Fragment
```
Reading DNA sequence from: ATCGAATTCGGGATCCAAA
Topology: circular
Single-cut behavior: intact circle (yields 1 fragment)
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

================================================================================
DIGESTION RESULTS
================================================================================
Mode: circular
Sequence length: 19 bp
Total cuts: 2
Cut positions: 4, 11
Fragments generated: 2

Fragment Details:
--------------------------------------------------------------------------------
Index    Start    End      Length     Wraps    Boundaries
--------------------------------------------------------------------------------
0        4        11       7          No       4(EcoRI) -> 11(BamHI)
1        11       4        12         Yes      11(BamHI) -> 4(EcoRI)
--------------------------------------------------------------------------------
✓ Fragment lengths sum correctly to 19 bp

Cut Site Details:
--------------------------------------------------------------------------------
  Position 4: EcoRI (site: GAATTC, overhang: 5' overhang)
  Position 11: BamHI (site: GGATCC, overhang: 5' overhang)
```

### Circular DNA - Single Cut (Default: Intact Circle)
```
Topology: circular
Single-cut behavior: intact circle (yields 1 fragment)
DNA sequence length: 100 bp

Fragments generated: 1

Fragment Details:
--------------------------------------------------------------------------------
Index    Start    End      Length     Wraps    Boundaries
--------------------------------------------------------------------------------
0        0        0        100        Yes      START -> END
--------------------------------------------------------------------------------
✓ Fragment lengths sum correctly to 100 bp
```

### Circular DNA - Single Cut (Linearization Mode)
```
Topology: circular
Single-cut behavior: linearizes (yields 2 fragments)
DNA sequence length: 100 bp

Fragments generated: 2

Fragment Details:
--------------------------------------------------------------------------------
Index    Start    End      Length     Wraps    Boundaries
--------------------------------------------------------------------------------
0        30       100      70         No       30(EcoRI) -> 30(EcoRI)
1        0        30       30         No       30(EcoRI) -> 30(EcoRI)
--------------------------------------------------------------------------------
✓ Fragment lengths sum correctly to 100 bp
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
├── sim.py                          # Main simulator script with linear/circular support
├── fragment_calculator.py          # Fragment computation module (linear, circular, restriction maps, gel simulation)
├── gel_ladders.py                  # DNA ladder presets for gel simulation
├── ligation_compatibility.py       # Ligation compatibility analysis (NEW!)
├── graphics.py                     # SVG/PNG graphics generation
├── enzymes.json                    # Extended enzyme database (350+ enzymes)
├── sample_dna.fasta                # Sample DNA sequence for testing
├── synthetic_restriction_test.fasta # Test sequence for synthetic enzymes
├── tests/
│   ├── test_circular.py                # Comprehensive circular DNA tests
│   ├── test_multi_enzyme.py            # Test cases for multi-enzyme functionality
│   ├── test_enhanced_features.py       # Test cases for IUPAC and enhanced features
│   ├── test_iupac.py                   # IUPAC expansion tests
│   ├── test_restriction_map.py         # Restriction map visualization tests
│   ├── test_gel_ascii.py               # ASCII gel simulation tests
│   ├── test_ligation_compatibility.py  # Ligation compatibility tests (NEW!)
│   └── test_ligation_sequence.fasta    # Test sequence for ligation analysis
├── prompts/
│   ├── prompt.txt                  # Original requirements
│   ├── prompt 2.txt                # Enhancement requirements
│   ├── prompt 3.txt                # Multi-enzyme requirements
│   ├── prompt 4.txt                # Advanced features requirements
│   ├── prompt 5                    # Implementation requirements
│   ├── prompt 6.txt                # Additional requirements
│   ├── prompt 7.txt                # Further requirements
│   ├── prompt 8.txt                # Circular DNA requirements
│   ├── prompt 9.txt                # Restriction map requirements
│   ├── prompt 10.txt               # Gel simulation requirements
│   └── prompt 11.txt               # Graphics output requirements
│   └── prompt 12.txt               # Ligation compatibility requirements (THIS IMPLEMENTATION)
└── README.md                       # This file
```

## Testing

### Manual Testing

Test the simulator with the included sample files:

```bash
# Linear DNA - single enzyme
python sim.py --seq sample_dna.fasta --enz EcoRI

# Linear DNA - multiple enzymes
python sim.py --seq sample_dna.fasta --enz EcoRI BamHI

# Circular DNA - two cuts (shows wrap-around)
python sim.py --seq sample_dna.fasta --enz EcoRI BamHI --circular

# Circular DNA - single cut (intact circle)
python sim.py --seq sample_dna.fasta --enz EcoRI --circular

# Circular DNA - single cut (linearized)
python sim.py --seq sample_dna.fasta --enz EcoRI --circular --circular_single_cut_linearizes

# Case-insensitive enzyme names
python sim.py --seq "ATCGAATTCGGGATCCAAA" --enz ecori bamhi

# Direct sequence with multiple enzymes
python sim.py --seq "ATCGAATTCGGGATCCAAA" --enz EcoRI BamHI HindIII
```

### Automated Testing

Run the comprehensive test suite with pytest:

```bash
# Run all tests
python -m pytest tests/ -v

# Run only circular DNA tests
python -m pytest tests/test_circular.py -v

# Run specific test class
python -m pytest tests/test_circular.py::TestCircularTwoCuts -v
```

**Test Coverage:**
- 26 circular DNA tests covering all edge cases
- Linear mode regression tests
- Boundary annotation tests
- Fragment validation tests
- Edge cases: zero-length sequences, modulo normalization, etc.

## Circular DNA Digestion

### Overview

Circular DNA (such as plasmids) requires special handling for fragment calculation because the DNA molecule has no "ends" - position 0 is adjacent to position L-1 (where L is the sequence length). This creates the possibility of **wrap-around fragments** that span the origin.

### Circular Mode Behavior

#### No Cuts
- Returns one intact circular molecule
- Length = full sequence length
- Marked as wrapping (conceptually)

#### One Cut
**Default behavior (`--circular` only):**
- Returns one fragment of length L
- Treats the molecule as intact (a single nick doesn't fragment a circle in silico unless specified)

**Linearization mode (`--circular --circular_single_cut_linearizes`):**
- Returns two fragments
- Effectively "opens" the circle at the cut site
- Fragment lengths sum to L

#### Two or More Cuts
- Returns N fragments (where N = number of cuts)
- N-1 fragments are regular intervals between consecutive cuts
- 1 fragment is a **wrap-around fragment** that spans from the last cut, around position 0, to the first cut
- The wrap-around fragment has:
  - `start > end` (e.g., start=90, end=10 in a 100bp sequence)
  - `wraps: true` flag
  - `length = (L - start) + end`

### Example: Circular Plasmid with Two Cuts

Given a 120 bp circular plasmid with cuts at positions 30 and 90:

**Linear mode would produce:**
- Fragment 1: [0-30] = 30 bp
- Fragment 2: [30-90] = 60 bp  
- Fragment 3: [90-120] = 30 bp
- Total: 3 fragments

**Circular mode produces:**
- Fragment 1: [30-90] = 60 bp (regular)
- Fragment 2: [90-30] = 60 bp (wrap-around: wraps=true)
  - Length = (120-90) + 30 = 60 bp
- Total: 2 fragments

### Boundary Annotations

Each fragment includes boundary information showing which enzyme(s) created each cut:

```
Boundaries:
  left_cut: {pos: 30, enzymes: [EcoRI]}
  right_cut: {pos: 90, enzymes: [BamHI]}
```

If multiple enzymes cut at the same position, all are listed in the boundary annotation.

## Technical Details

### How It Works

1. **Database Loading**: Loads enzyme database from `enzymes.json` if present, otherwise uses built-in database
2. **Sequence Input**: Reads DNA sequence from string or file, converting to uppercase
3. **Topology Selection**: Determines linear or circular mode based on `--circular` flag
4. **Enzyme Validation**: Validates and normalizes enzyme names (case-insensitive) with similarity matching for errors
5. **IUPAC Expansion**: Converts IUPAC degenerate bases to regex character classes for pattern matching
6. **Pattern Matching**: Uses lookahead regex to find all overlapping recognition sites
7. **Cut Position Calculation**: Computes break positions using cut_index specifications
8. **Metadata Collection**: Collects enzyme name, site, cut_index, and overhang_type for each cut
9. **Cut Position Merging**: Combines cut positions from all enzymes, deduplicating overlapping cuts
10. **Fragment Calculation**: 
    - Linear mode: Standard interval calculation
    - Circular mode: Includes wrap-around fragment calculation
11. **Boundary Annotation**: Maps enzymes to fragment boundaries based on cut positions
12. **Validation**: Verifies fragment lengths sum to sequence length
13. **Output Generation**: Displays comprehensive results including positions, wrapping status, and boundary enzymes

### Code Structure

**Main Module (sim.py):**
- `iupac_to_regex()`: Converts IUPAC degenerate bases to regex character classes
- `normalize()`: Normalizes enzyme names by removing diacritics and converting to lowercase
- `load_enzyme_database()`: Loads enzyme information from JSON file with overhang_type support
- `read_dna_sequence()`: Handles input from string or file with validation
- `find_cut_sites()`: Locates all recognition sequences using IUPAC expansion
- `find_cut_positions_linear()`: Finds cut positions for a single enzyme
- `merge_cut_positions()`: Combines and deduplicates cut positions from multiple enzymes
- `find_closest_enzyme_names()`: Provides helpful suggestions for unknown enzyme names
- `main()`: Orchestrates the entire process with linear/circular mode support

**Fragment Calculator Module (fragment_calculator.py):**
- `compute_fragments()`: Core algorithm for fragment calculation supporting both topologies
  - Handles 0, 1, or N cuts
  - Implements wrap-around logic for circular DNA
  - Attaches boundary annotations with enzyme metadata
- `validate_fragment_total()`: Ensures fragment lengths sum to sequence length
- `build_restriction_map()`: Generates text-mode restriction maps
- `simulate_gel()`: Generates ASCII agarose gel electrophoresis visualization
- `gel_coefficients()`, `calculate_migration_row()`: Migration model for gel simulation
- `merge_bands()`, `get_band_glyph()`: Band rendering and merging logic
- `render_circular_topology_bands()`: Special rendering for SC/OC plasmid forms

**Gel Ladders Module (gel_ladders.py):**
- `get_ladder()`: Returns fragment sizes for standard DNA ladders
- `get_available_ladders()`: Lists available ladder presets
- Presets: 100bp, 1kb, broad

**Ligation Compatibility Module (ligation_compatibility.py):**
- `revcomp()`: Calculate reverse complement of DNA sequences
- `calculate_gc_percent()`: Calculate GC percentage of overhang sequences
- `calculate_tm()`: Estimate melting temperature using Wallace rule
- `are_compatible()`: Check if two fragment ends are compatible for ligation
- `is_directional()`: Determine if a pair is directional (non-palindromic)
- `calculate_compatibility()`: Analyze all fragment ends for compatibility
- `format_pairs_output()`: Format results as list of compatible pairs
- `format_matrix_output()`: Format results as compatibility matrix
- `format_detailed_output()`: Format results with full details and heuristics
- `export_to_json()`: Export compatibility results to JSON
- `analyze_enzyme_pair_theoretical()`: Analyze theoretical enzyme compatibility

## Future Enhancements

Possible future enhancements may include:
- Export to standard formats (GenBank, CSV)
- Support for more complex enzyme behaviors (e.g., methylation sensitivity, temperature requirements)
- Enhanced Type IIS enzyme support with better visualization
- Dam/Dcm methylation blocking
- Star activity prediction
- Theoretical enzyme pair compatibility (query without digest)
- Multi-step cloning strategy planning

**Recently Implemented:**
- ✅ Text-mode restriction map visualization (prompt 9)
- ✅ ASCII agarose gel simulation (prompt 10)
- ✅ Ligation compatibility analysis (prompt 12)

## License

This project is open source and available for educational use.

## Contributing

This is a learning project designed for beginners. Feel free to use and modify for educational purposes!
