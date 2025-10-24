# Fragment Sequence Display Feature

## Overview

The Fragment Sequence Display feature allows users to view the actual DNA sequences for each restriction enzyme digest fragment. This enhancement provides a complete view of fragment information including:

- Fragment positions and length
- DNA sequence (5'→3' orientation)
- Enzyme cut site information
- Overhang types and sequences
- End information for ligation planning

## Implementation

### iOS App (Swift)

#### Files Modified/Created

1. **FragmentList.swift** - Enhanced fragment list view
   - Now displays sequence previews in each fragment row
   - Shows overhang information for both ends
   - Links to detailed fragment view
   - Accepts `fullSequence` parameter (for future enhancements)

2. **FragmentDetailView.swift** - New detailed fragment view
   - Complete fragment information card
   - Detailed end information (5' and 3')
   - Full sequence display with formatting
   - Expandable view for long sequences
   - Copy-to-clipboard functionality
   - Position markers for sequence navigation

3. **HomeView.swift** - Updated digest engine call
   - Changed `returnSequences: false` to `returnSequences: true`
   - Passes full sequence to FragmentList
   - Enables sequence retrieval in digest results

#### Features

##### Fragment List View
- **Sequence Preview**: Shows first 20 and last 17 bases for fragments > 40bp
- **Overhang Display**: Shows enzyme names and sticky sequences for both ends
- **Navigation**: Tap any fragment to view full details

##### Fragment Detail View
- **Summary Card**: Fragment length and position information
- **Ends Information**: 
  - Enzyme names
  - Overhang types (5', 3', or blunt)
  - Sticky sequences for ligation planning
- **Sequence Display**:
  - Formatted with position markers
  - Grouped in 10-base increments
  - Show/Hide toggle for long sequences (>120bp)
  - Selectable text for easy copying
  - Dedicated copy button

### Python Reference Implementation

#### File: `scripts/fragment_sequences.py`

Provides reference implementation demonstrating:

1. **Core Function**: `extract_fragment_sequence()`
   - Extracts DNA sequence between two positions
   - Handles circular DNA with wrap-around
   - Returns 5'→3' orientation

2. **Complete Workflow**: `get_fragment_with_sequence()`
   - Performs restriction digest
   - Returns fragments with full information
   - Includes end metadata and sequences

3. **Formatting**: `format_sequence_display()`
   - Pretty-prints DNA sequences
   - Position markers and grouping
   - Configurable line length and group size

4. **Examples**: Three working examples
   - Linear DNA with EcoRI
   - Circular plasmid with multiple enzymes
   - Test sequence from data files

#### Usage Example

```python
from fragment_sequences import get_fragment_with_sequence
import json

# Load enzyme database
with open('data/enzymes.json', 'r') as f:
    enzyme_list = json.load(f)
    enzyme_db = {e['name']: e for e in enzyme_list}

# Perform digest
dna = "ATGCGAATTCGCTAGCTAGC"
enzymes = ["EcoRI"]
fragments = get_fragment_with_sequence(dna, enzymes, enzyme_db, circular=False)

# Access fragment information
for frag in fragments:
    print(f"Fragment {frag['index'] + 1}:")
    print(f"  Position: [{frag['start']} .. {frag['end']})")
    print(f"  Length: {frag['length']} bp")
    print(f"  Sequence: {frag['sequence']}")
    if frag['left_end']:
        print(f"  5' End: {frag['left_end']['enzyme']} ({frag['left_end']['overhang_type']})")
    if frag['right_end']:
        print(f"  3' End: {frag['right_end']['enzyme']} ({frag['right_end']['overhang_type']})")
```

## User Experience

### Viewing Fragment Sequences

1. **Enter DNA sequence** in the main view
2. **Select enzymes** for digestion
3. **Tap "Digest"** to perform restriction analysis
4. **Tap "View Fragments"** to see all fragments
5. **Tap any fragment** to view detailed information including full sequence

### Fragment List

Each fragment row displays:
- Fragment length (e.g., "101 bp")
- Position range (e.g., "[0 .. 101)")
- Enzyme information for 5' and 3' ends
- Sequence preview with truncation indicator

### Fragment Detail View

The detail view provides:
- **Fragment Information Card**
  - Total length in base pairs
  - Position in original sequence
  
- **Fragment Ends Card**
  - 5' End (Left): Enzyme, overhang type, sticky sequence
  - 3' End (Right): Enzyme, overhang type, sticky sequence
  
- **DNA Sequence Card**
  - Position-numbered sequence display
  - 10-base grouping for readability
  - Show All/Show Less toggle for long sequences
  - Copy to clipboard button

## Technical Details

### Sequence Retrieval

The Swift `DigestEngine` class supports optional sequence retrieval via the `returnSequences` parameter:

```swift
let options = DigestOptions(
    circular: circular, 
    returnSequences: true  // Enable sequence extraction
)
let fragments = engine.digest(options: options)
```

When enabled, each `Fragment` object includes:
- `sequence: String?` - The DNA sequence in 5'→3' orientation
- Sequence slicing handles both linear and circular topologies
- Circular fragments may wrap around the origin

### Sequence Formatting

Sequences are formatted for readability:
- **Position Markers**: Every line starts with base position (1-indexed)
- **Grouping**: Bases grouped in tens (e.g., "ATGCGAATTC GCTAGCTAG")
- **Truncation**: Long sequences show first/last portions with "..." indicator
- **Monospaced Font**: Uses system monospaced design for alignment

### Memory Considerations

For large sequences:
- Sequences only retrieved when `returnSequences: true`
- Fragment list shows truncated previews (40bp max)
- Detail view allows toggling full sequence display
- Selectable text enables external analysis without copying entire sequence

## Integration with Existing Features

### Export Functionality
- CSV export includes fragment positions (sequences can be extracted)
- GenBank export includes fragment features with positions
- Future enhancement: Add sequence column to CSV export

### Gel View
- Fragments maintain sequence data throughout analysis
- Future enhancement: Show sequence in gel band tooltips

### Map View (Circular DNA)
- Fragment positions already visualized
- Future enhancement: Sequence display on fragment hover

## Future Enhancements

1. **Sequence Search**: Find specific sequences within fragments
2. **Translation**: Show amino acid translation for coding sequences
3. **GC Content**: Calculate and display GC percentage per fragment
4. **Restriction Sites**: Identify additional restriction sites within fragments
5. **Export Sequences**: FASTA export for individual or all fragments
6. **Reverse Complement**: Toggle to show reverse complement strand
7. **Highlight Features**: Mark specific motifs or patterns in sequence

## Testing

### Test Cases

1. **Linear DNA, Single Enzyme**
   - Verify fragment count matches expected cuts + 1
   - Check sequence boundaries match positions
   - Confirm overhang information is correct

2. **Circular DNA, Multiple Enzymes**
   - Test wrap-around fragment sequences
   - Verify circular topology handling
   - Check enzyme ordering at boundaries

3. **Long Sequences**
   - Test truncation in list view
   - Verify Show All/Show Less toggle
   - Check copy functionality

4. **Edge Cases**
   - No cuts (single fragment = full sequence)
   - Single cut on circular DNA
   - Blunt end cutters (no sticky sequence)
   - Multiple enzymes cutting at same position

### Python Reference Tests

Run the reference implementation:

```bash
cd scripts
python3 fragment_sequences.py
```

Expected output: Three examples with detailed fragment information.

## Accessibility

The feature includes full accessibility support:

- **VoiceOver Labels**: All fragments have descriptive labels
- **Accessibility Hints**: Guidance for navigation and actions
- **Dynamic Type**: Text scales with system font size settings
- **Selectable Text**: Sequences can be selected by assistive technologies
- **Semantic Structure**: Proper heading hierarchy and grouping

## Performance

- **Lazy Loading**: Sequences only computed when requested
- **Efficient Slicing**: O(n) sequence extraction where n = fragment length
- **View Optimization**: List uses NavigationLink for on-demand detail loading
- **Memory**: Sequences stored in Fragment objects, released when view dismissed

## Dependencies

- **DigestCore**: Swift package providing digest engine
- **SwiftUI**: UI framework for views
- **Foundation**: String and data manipulation

Python reference requires:
- `fragment_calculator.py`: Fragment computation
- `sim.py`: Enzyme cut site finding
- `data/enzymes.json`: Enzyme database

## Conclusion

The Fragment Sequence Display feature provides users with complete information about their restriction digest results. The implementation leverages existing digest engine capabilities and presents data in an intuitive, accessible interface. The Python reference implementation ensures consistency and provides a foundation for future enhancements.

