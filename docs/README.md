# Genomancer - DNA Restriction Analysis for iOS

<p align="center">
  <img src="../../Screenshots/Intro.png" alt="Genomancer App" width="250"/>
</p>

**Genomancer** is a powerful, elegant iOS app for simulating restriction enzyme digestion of DNA sequences. Whether you're a molecular biology researcher, student, or educator, Genomancer provides professional-grade analysis tools in an intuitive, mobile-friendly interface.

[![Platform](https://img.shields.io/badge/platform-iOS-lightgrey.svg)](https://www.apple.com/ios/)
[![Swift](https://img.shields.io/badge/Swift-5.9-orange.svg)](https://swift.org)
[![License](https://img.shields.io/badge/license-Open%20Source-blue.svg)](../../LICENSE)

## Features

<<<<<<< HEAD
### 🧬 Comprehensive Restriction Enzyme Analysis
- **350+ Enzyme Database**: Complete library of restriction enzymes with recognition sites, cut positions, and overhang types
- **Linear & Circular DNA Support**: Analyze both linear DNA fragments and circular plasmids
- **Multiple Enzyme Digests**: Select and combine multiple enzymes for complex digest simulations
- **IUPAC Support**: Full support for degenerate base codes (N, R, Y, S, W, K, M, B, D, H, V)

### 📊 Visualization Tools
- **Fragment List**: Detailed view of all DNA fragments with sizes, positions, and end characteristics
- **Gel Electrophoresis Simulation**: Realistic agarose gel visualization showing fragment migration
- **Plasmid Map**: Beautiful circular maps for plasmid digests with labeled cut sites
- **Ligation Compatibility Analysis**: Analyze which fragment ends can be ligated together

### 🎯 Powerful Analysis Features
- **Fragment Sequences**: View exact DNA sequences for each fragment with overhang information
- **Overhang Analysis**: Detailed information about 5' overhangs, 3' overhangs, and blunt ends
- **Cut Site Details**: Complete information about recognition sites and cut positions
- **Wrap-Around Fragments**: Correct handling of circular DNA fragments that span the origin

### 📁 Flexible Input Options
- **FASTA Import**: Import DNA sequences directly from FASTA files
- **Manual Entry**: Type or paste sequences directly into the app
- **Format Detection**: Automatic detection and parsing of FASTA format
- **Sequence Validation**: Real-time validation of DNA sequence characters

## Screenshots

<p align="center">
  <img src="../../Screenshots/Digest Intro.png" alt="Digest View" width="200"/>
  <img src="../../Screenshots/Enzymes Intro.png" alt="Enzyme Selection" width="200"/>
  <img src="../../Screenshots/Fragments.png" alt="Fragment List" width="200"/>
  <img src="../../Screenshots/Gel.png" alt="Gel Simulation" width="200"/>
</p>

<p align="center">
  <img src="../../Screenshots/Map.png" alt="Plasmid Map" width="200"/>
  <img src="../../Screenshots/Ligation Analysis.png" alt="Ligation Analysis" width="200"/>
</p>

## How to Use

### Basic Workflow

1. **Enter Your DNA Sequence**
   - Type or paste your sequence directly
   - Import a FASTA file using the "Import FASTA File" button
   - Supports both raw sequences and FASTA format with headers

2. **Select Topology**
   - Toggle "Circular (plasmid)" ON for circular DNA (plasmids, viral genomes)
   - Leave OFF for linear DNA (PCR products, genomic fragments)

3. **Choose Restriction Enzymes**
   - Tap "Choose Enzyme(s)" to browse the enzyme database
   - Search for specific enzymes by name
   - Select multiple enzymes for double or triple digests
   - View selected enzymes as chips with easy removal

4. **Run the Digest**
   - Tap the "Digest" button to perform the analysis
   - Wait for results to be calculated

5. **Explore Results**
   - **View Fragments**: See all generated fragments with detailed information
   - **View Map**: Circular plasmid map showing cut sites (circular DNA only)
   - **View Gel**: Simulated agarose gel electrophoresis
   - **Ligation Analysis**: Check which fragment ends are compatible for ligation
=======
- **Linear and Circular DNA Support**: Full support for both linear and circular DNA topology with wrap-around fragment calculation
- **Fragment Sequence Extraction**: Returns actual DNA sequences for each fragment with detailed overhang analysis.
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

**Core functionality requires no dependencies!** This simulator uses only Python standard library modules:
- `argparse` for command-line arguments
- `json` for loading enzyme database
- `re` for pattern matching and IUPAC expansion
- `sys` for system operations
- `unicodedata` for name normalization
- `typing` for type hints

All main features (digest, analysis, SVG output) work out of the box with Python 3.7+.

### Optional Dependencies

For **PNG export** from SVG graphics (when using `--png` flag):
```bash
pip install cairosvg
```

Note: PNG export requires both `cairo` system library and `cairosvg` Python package. SVG output works without any installation.

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

**Field validation:**
- `cut_index` must satisfy: `0 ≤ cut_index ≤ len(site)` (cuts within or at boundaries of recognition site)
- `overhang_type` must be one of: `"5' overhang"`, `"3' overhang"`, `"Blunt"`, or `"Unknown"`
- `site` must contain only valid IUPAC nucleotide codes (case-insensitive)

**Note:** The simulator validates `enzymes.json` on load and will skip invalid entries with a warning. Overhang type is displayed verbatim (normalized) in output.

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

##### Duplicate enzyme handling:

When you specify the same enzyme multiple times, the simulator automatically assigns unique display labels with `#2`, `#3` suffixes for clarity:

```bash
# Specify EcoRI twice to track different recognition sites separately
python sim.py --seq "GAATTCGGGGGAATTC" --enz EcoRI EcoRI --print-map
```

**Output will show:**
```
Enzyme: EcoRI
Site:   GAATTC
...

Enzyme: EcoRI#2
Site:   GAATTC
...

Cut Site Details:
  Position 1: EcoRI ...
  Position 12: EcoRI#2 ...
```

This is particularly useful when:
- Tracking multiple instances of the same enzyme in restriction maps
- Distinguishing between different cut sites in the same sequence
- Creating complex multi-enzyme digestion strategies

**Note:** The actual enzyme behavior is identical (same recognition site and cut pattern), but the display labels help you identify which occurrence you're looking at in the output.

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

### Overhang Semantics and Calculation

The simulator uses a **centralized overhang calculation function** (`compute_end_metadata`) to ensure consistency across all outputs (console, CSV, GenBank, JSON). This function is the single source of truth for overhang metadata.

#### How Overhang Length is Calculated

For standard Type II restriction enzymes with palindromic recognition sites, the overhang length is calculated as:

```
overhang_length = |2 × cut_index - site_length|
```

**Examples:**
- **EcoRI** (GAATTC, cut_index=1): |2×1 - 6| = |2-6| = **4 bp** ✓
- **HindIII** (AAGCTT, cut_index=1): |2×1 - 6| = **4 bp** ✓
- **PstI** (CTGCAG, cut_index=5): |2×5 - 6| = |10-6| = **4 bp** ✓
- **SmaI** (CCCGGG, cut_index=3): |2×3 - 6| = |6-6| = **0 bp** (blunt) ✓

#### Overhang Orientation and Reporting

All sticky end sequences are reported in **5'→3' orientation** relative to the fragment strand:

#### 5' Overhang (EcoRI)
```
Recognition site: GAATTC
Cut position: G^AATTC (cuts after position 1)
Overhang type: 5' overhang
Overhang length: 4 bp

Top strand:    5'---G AATT---3'
Bottom strand: 3'---CTTAA G---5'
                      ^^^^
                    4 bp 5' overhang

Reported sticky sequence: AATT (5'→3')
```

**In output:**
- Console: "5' overhang, 4 bp: AATT"
- CSV: `left_overhang_len=4, left_end_bases=AATT`
- GenBank: `/note="overhang=5' overhang; k=4"`

#### 3' Overhang (PstI)
```
Recognition site: CTGCAG
Cut position: CTGCA^G (cuts after position 5)
Overhang type: 3' overhang
Overhang length: 4 bp

Top strand:    5'---CTGCA G---3'
Bottom strand: 3'---G ACGTC---5'
                    ^^^^
                  4 bp 3' overhang

Reported sticky sequence: TGCA (5'→3' on the protruding strand)
```

**In output:**
- Console: "3' overhang, 4 bp: TGCA"
- CSV: `right_overhang_len=4, right_end_bases=TGCA`
- GenBank: `/note="overhang=3' overhang; k=4"`

#### Blunt Ends (SmaI)
```
Recognition site: CCCGGG
Cut position: CCC^GGG (cuts at position 3, middle)
Overhang type: Blunt
Overhang length: 0 bp

Top strand:    5'---CCC GGG---3'
Bottom strand: 3'---GGG CCC---5'
                    ^
                  No overhang

Reported sticky sequence: "" (empty)
```

**In output:**
- Console: "Blunt, 0 bp"
- CSV: `overhang_len=0, end_bases=""`
- GenBank: `/note="overhang=Blunt; k=0"`

#### Output Consistency Guarantee

The centralized `compute_end_metadata()` function ensures that:
- Overhang lengths are **always** calculated correctly (EcoRI = 4 bp, not 5 or 0)
- Sticky sequences are **always** reported in the same orientation
- Console, CSV, GenBank, and JSON outputs **always** match
- Left and right fragment ends use the **same** calculation logic

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

The `--lanes-config` option allows you to define multiple lanes with different enzyme combinations in a single gel simulation. This is useful for comparing different digestion strategies side-by-side.

**Option:**
- `--lanes-config <json>`: Path to JSON file or inline JSON string defining multiple lanes

**Lane Configuration Format:**

Each lane in the configuration is a JSON object with the following optional fields:

```json
{
  "label": "Lane Label",           // Display name for the lane
  "enzymes": ["EcoRI", "HindIII"], // List of enzymes (optional, uses global --enz if omitted)
  "circular": true,                // Override global topology for this lane (optional)
  "notes": "Optional description"  // Additional notes (optional)
}
```

**Example lanes.json file:**

```json
[
  {"label": "Uncut", "enzymes": [], "circular": true},
  {"label": "EcoRI", "enzymes": ["EcoRI"], "circular": true},
  {"label": "EcoRI+HindIII", "enzymes": ["EcoRI", "HindIII"], "circular": true}
]
```

**Usage:**

```bash
# Using a JSON file
python sim.py --seq plasmid.fasta --lanes-config lanes.json --gel-only

# Using inline JSON (single lane example)
python sim.py --seq plasmid.fasta --lanes-config '[{"label":"Test","enzymes":["EcoRI"]}]' --gel-only

# Multiple lanes comparing different enzyme combinations
cat > lanes.json <<'JSON'
[
  {"label":"Uncut","enzymes":[],"circular":true},
  {"label":"EcoRI","enzymes":["EcoRI"],"circular":true},
  {"label":"EcoRI+HindIII","enzymes":["EcoRI","HindIII"],"circular":true}
]
JSON

python sim.py --seq "ATCGAATTCGGGATCCAAA" --lanes-config lanes.json --gel-only
```

**Notes:**
- If `enzymes` is omitted from a lane configuration, it will use the enzymes specified via `--enz`
- The `circular` field in a lane config overrides the global `--circular` flag for that specific lane
- Each lane is processed independently and fragments are calculated for each enzyme combination

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
>>>>>>> 351894ba3fa48a16c78d978c1f96c24f2012fa30

### Example Use Cases

#### Verifying a Plasmid Digest
1. Import your plasmid sequence (e.g., pUC19)
2. Enable "Circular (plasmid)" mode
3. Select restriction enzymes (e.g., EcoRI, BamHI)
4. Run digest and check fragment sizes match expected results

#### Planning a Cloning Strategy
1. Enter your insert and vector sequences
2. Test different enzyme combinations
3. Use Ligation Analysis to verify compatible ends
4. Check the gel simulation to predict band patterns

#### Teaching Molecular Biology
1. Import example sequences
2. Demonstrate restriction digests with different enzymes
3. Show students gel migration patterns
4. Explain overhang types and ligation compatibility

## Detailed Features

### Fragment Analysis
Each fragment includes:
- **Position**: Start and end coordinates (0-based indexing)
- **Length**: Fragment size in base pairs
- **Sequence**: Complete DNA sequence (5' to 3')
- **Boundaries**: Which enzymes created each end
- **Overhang Type**: 5' overhang, 3' overhang, or blunt
- **Overhang Sequence**: Exact bases in sticky ends

### Ligation Compatibility
The ligation analysis tool checks:
- **Overhang Compatibility**: Whether two ends have matching sticky sequences
- **Overhang Type Matching**: Both ends must be 5' or both 3'
- **Directionality**: Identifies non-palindromic overhangs for directional cloning
- **GC Content**: Calculates GC% of overhang sequences
- **Melting Temperature**: Estimates Tm using Wallace rule

### Gel Simulation
Realistic gel electrophoresis simulation featuring:
- **Size-Based Migration**: Fragments migrate based on molecular weight
- **Agarose Percentage**: Optimized for different size ranges
- **DNA Ladder**: Standard marker for size comparison
- **Band Intensity**: Visual representation of DNA concentration
- **Circular DNA Forms**: Shows supercoiled (SC) and open circular (OC) forms

### Plasmid Map
Beautiful circular maps showing:
- **Cut Site Markers**: Radial lines indicating restriction sites
- **Position Labels**: Base pair positions for each cut
- **Enzyme Names**: Clear labeling of which enzyme cuts where
- **Tick Marks**: Scale indicators every 1000 bp
- **Multiple Cuts**: Visual indication when multiple enzymes cut at same position

## Enzyme Database

Genomancer includes a comprehensive database of **350+ restriction enzymes**, including:

### Common Enzymes
- **EcoRI** (G^AATTC) - 5' overhang, 4 bp
- **BamHI** (G^GATCC) - 5' overhang, 4 bp
- **HindIII** (A^AGCTT) - 5' overhang, 4 bp
- **PstI** (CTGCA^G) - 3' overhang, 4 bp
- **NotI** (GC^GGCCGC) - 5' overhang, 4 bp
- **EcoRV** (GAT^ATC) - Blunt end

### Enzyme Categories
- **Type II Enzymes**: Standard restriction enzymes with palindromic recognition sites
- **Type IIS Enzymes**: Cut outside their recognition sequence (e.g., BsaI, BsmBI)
- **Blunt Cutters**: Produce blunt ends (e.g., SmaI, EcoRV)
- **4-Base Cutters**: Frequent cutters (e.g., TaqI, MseI)
- **6-Base Cutters**: Standard molecular biology enzymes
- **8-Base Cutters**: Rare cutters (e.g., NotI, AscI)

Each enzyme entry includes:
- Recognition sequence (with IUPAC support)
- Cut position index
- Overhang type (5', 3', or blunt)
- Common usage notes

## Technical Details

### System Requirements
- **iOS 16.0 or later**
- Compatible with iPhone and iPad
- Optimized for all screen sizes
- Supports Dark Mode
- VoiceOver accessible

### Architecture
- Built with **SwiftUI** for modern, responsive UI
- **DigestCore** Swift package for enzyme calculations
- Native performance with no external dependencies
- Efficient handling of large sequences (tested up to 50 kb)

### Sequence Handling
- **Maximum Length**: Tested with sequences up to 50,000 bp
- **Format Support**: Raw sequences and FASTA format
- **Character Validation**: Real-time validation of IUPAC codes
- **Memory Efficient**: Lazy loading and optimized fragment storage

### Enzyme Matching Algorithm
- **IUPAC Expansion**: Full support for degenerate bases
- **Overlapping Sites**: Correctly handles overlapping recognition sequences
- **Cut Position Calculation**: Precise positioning including Type IIS enzymes
- **Circular Wrap-Around**: Correct fragment calculation across origin

## Privacy & Data

**Genomancer respects your privacy:**
- ✅ No data collection
- ✅ No analytics or tracking
- ✅ No internet connection required
- ✅ All processing happens locally on your device
- ✅ Your sequences never leave your iPhone/iPad

## Support & Feedback

### Having Issues?
- Check that your sequence contains only valid IUPAC DNA codes (A, T, C, G, N, R, Y, S, W, K, M, B, D, H, V)
- Ensure you've selected at least one enzyme before digesting
- For circular DNA features, make sure "Circular (plasmid)" toggle is enabled
- Try reimporting FASTA files if sequence isn't recognized

### Contact
- **Email**: flavorislab@gmail.com

### Contributing
We welcome contributions:
- Bug fixes and improvements
- New enzyme additions
- Feature requests and ideas
- Documentation improvements

## FAQ

### Q: Can I use my own enzyme database?
**A:** Currently, the app uses the built-in enzyme database.

### Q: What's the difference between circular and linear mode?
**A:** Circular mode treats the DNA as a closed loop (like a plasmid), which affects fragment calculation. In circular DNA with 2+ cuts, the last fragment wraps around from the end back to the beginning. Linear mode treats DNA as having distinct ends.

### Q: Why do I see "wrap-around" fragments?
**A:** In circular DNA (plasmids), if there are multiple restriction sites, one fragment will span from the last cut site, around position 0, back to the first cut site. This is biologically accurate for circular molecules.

### Q: Can I export my results?
**A:** Results can be viewed and analyzed within the app. Export functionality (CSV, GenBank, images) is planned for future updates.

### Q: Does this replace lab work?
**A:** No! Genomancer is a planning and verification tool. Always confirm computational predictions with actual lab experiments. Enzyme activity can be affected by temperature, buffer conditions, methylation, and other factors not modeled by the app.

### Q: What are the overhang numbers (5', 3', blunt)?
**A:** These indicate the type of DNA ends created:
- **5' overhang**: Top strand extends beyond bottom strand (e.g., EcoRI)
- **3' overhang**: Bottom strand extends beyond top strand (e.g., PstI)
- **Blunt**: Both strands cut at the same position (e.g., SmaI)

### Q: Why can't certain fragments ligate together?
**A:** For ligation, fragments must have:
- Compatible overhang types (both 5' OR both 3')
- Matching overhang lengths
- Complementary sticky-end sequences
- Or both must be blunt ends (with T4 DNA ligase)

## Related Resources

### Scientific Background
- [Restriction Enzymes (Wikipedia)](https://en.wikipedia.org/wiki/Restriction_enzyme)
- [REBASE - The Restriction Enzyme Database](http://rebase.neb.com/)
- [NEB's Restriction Enzyme Resource](https://www.neb.com/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities)

### Molecular Biology Tools
- [SnapGene](https://www.snapgene.com/) - Commercial DNA analysis software
- [Benchling](https://www.benchling.com/) - Cloud-based molecular biology platform
- [ApE](https://jorgensen.biology.utah.edu/wayned/ape/) - A plasmid Editor

### Educational Resources
- [Addgene's Plasmid Guide](https://www.addgene.org/protocols/plasmids-101/)
- [NCBI Molecular Biology Resources](https://www.ncbi.nlm.nih.gov/guide/molecular-biology/)

## Acknowledgments

Genomancer is built with:
- **Swift & SwiftUI** - Modern iOS development
- **DigestCore** - Custom Swift package for DNA analysis
- **NEB REBASE** - Enzyme database information

## Version History

### Version 1.0 (Current)
- Initial release
- 350+ enzyme database
- Linear and circular DNA support
- Fragment analysis
- Gel simulation
- Plasmid map visualization
- Ligation compatibility analysis
- FASTA file import
- Dark mode support
- Accessibility features

---

**Made with 🧬 for molecular biologists everywhere**