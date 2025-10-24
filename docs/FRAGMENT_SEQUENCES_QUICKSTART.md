# Fragment Sequences - Quick Start Guide

## What's New

The Genomancer app now displays the actual DNA sequence for each restriction enzyme fragment! 

## iOS App Usage

### Viewing Fragment Sequences

1. **Enter your DNA sequence** in the main view (FASTA or raw format)
2. **Select restriction enzymes** to use for digestion
3. **Tap "Digest"** to perform the analysis
4. **Tap "View Fragments"** to see all fragments

### In the Fragment List

Each fragment now shows:
- **Length** (e.g., "101 bp")
- **Position** range (e.g., "[0 .. 101)")
- **Enzyme cuts** at 5' and 3' ends with overhang sequences
- **Sequence preview** showing first and last bases

### Viewing Full Fragment Details

**Tap any fragment** to see:
- Complete fragment information
- Detailed end analysis (overhang types, sticky sequences)
- Full DNA sequence with position markers
- Copy to clipboard button

### Features in Detail View

- **Formatted Sequence**: 10-base grouping with line numbers
- **Show All/Show Less**: Toggle for long sequences (>120 bp)
- **Selectable Text**: Copy portions or use search
- **Copy Button**: One-tap sequence copy to clipboard

## Python Reference

A complete Python reference implementation is available in `scripts/fragment_sequences.py`.

### Quick Test

```bash
cd scripts
python3 fragment_sequences.py
```

This will run three examples:
1. Linear DNA with EcoRI
2. Circular plasmid with EcoRI and BamHI  
3. Real test sequence with HaeIII and HinfI

### Using in Your Scripts

```python
from fragment_sequences import get_fragment_with_sequence
import json

# Load enzymes
with open('data/enzymes.json', 'r') as f:
    enzyme_list = json.load(f)
    enzyme_db = {e['name']: e for e in enzyme_list}

# Digest DNA
dna = "ATGCGAATTCGCTAGC"
enzymes = ["EcoRI"]
fragments = get_fragment_with_sequence(dna, enzymes, enzyme_db, circular=False)

# Access sequences
for frag in fragments:
    print(f"Fragment {frag['index'] + 1}: {frag['sequence']}")
```

## Technical Details

### How It Works

The app uses the `DigestEngine` with `returnSequences: true` to extract the actual DNA sequence for each fragment based on:
- Start and end positions
- DNA topology (linear or circular)
- Wrap-around handling for circular DNA

### What's Included in Fragment Data

Each fragment contains:
- **Sequence**: The actual DNA bases (5'→3')
- **Position**: Start and end in original sequence
- **Length**: Fragment size in base pairs
- **Ends**: Information about both 5' and 3' ends
  - Source enzyme name
  - Overhang type (5', 3', or blunt)
  - Sticky sequence for ligation planning

## Examples

### Example 1: EcoRI Digest

```
Input: ATGCGAATTCGCTAGC (16 bp)
Enzyme: EcoRI (cuts G^AATTC)
Cut position: 5

Fragments:
1. ATGCG (5 bp)          - Natural 5' end, EcoRI 3' end (AATT overhang)
2. AATTCGCTAGC (11 bp)   - EcoRI 5' end (AATT overhang), Natural 3' end
```

### Example 2: Multiple Enzymes on Circular DNA

```
Input: 44 bp circular plasmid
Enzymes: EcoRI, BamHI
Cuts: EcoRI at 5, BamHI at 17, EcoRI at 36

Fragments:
1. 12 bp  [5..17)    EcoRI → BamHI
2. 19 bp  [17..36)   BamHI → EcoRI
3. 13 bp  [36..5)    EcoRI → EcoRI (wraps around)
```

## File Changes

### New Files
- `apps/ios/Genomancer/Genomancer/App/FragmentDetailView.swift` - Detail view with full sequence
- `scripts/fragment_sequences.py` - Python reference implementation
- `docs/FRAGMENT_SEQUENCES_FEATURE.md` - Complete feature documentation

### Modified Files
- `apps/ios/Genomancer/Genomancer/App/FragmentList.swift` - Added sequence previews
- `apps/ios/Genomancer/Genomancer/App/HomeView.swift` - Enabled sequence retrieval

## Tips

1. **Long Sequences**: Use "Show All" toggle to expand sequences longer than 120 bp
2. **Copy Specific Regions**: Select text directly from the formatted sequence view
3. **Ligation Planning**: Check sticky sequences to find compatible fragment ends
4. **Position Reference**: Use the line numbers to locate specific regions

## Future Enhancements

Potential additions:
- FASTA export for fragments
- GC content calculation
- Find/search within sequences
- Translation to amino acids
- Highlight restriction sites within fragments

## Need Help?

See the full documentation in `docs/FRAGMENT_SEQUENCES_FEATURE.md` for:
- Complete technical details
- API reference
- Testing guidelines
- Accessibility features
- Performance considerations

