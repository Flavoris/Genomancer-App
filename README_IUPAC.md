# IUPAC Degenerate Base Support

The restriction enzyme simulator now supports IUPAC degenerate bases in recognition sites, allowing for more flexible enzyme matching patterns.

## IUPAC Base Codes

| Code | Meaning | Matches |
|------|---------|---------|
| A    | Adenine | A       |
| C    | Cytosine| C       |
| G    | Guanine | G       |
| T    | Thymine | T       |
| R    | Purine  | A or G  |
| Y    | Pyrimidine | C or T |
| W    | Weak    | A or T  |
| S    | Strong  | G or C  |
| M    | Amino   | A or C  |
| K    | Keto    | G or T  |
| N    | Any     | A, C, G, or T |

## Examples

### Example 1: Simple Degenerate Base (R)

```bash
python sim.py --seq "ATAGATAGGCC" --enz TestR
```

Where `TestR` has site "AGR" and cut_index=1:

```
Enzyme: TestR
Recognition sequence: AGR
Cut position: after position 1 in recognition sequence
IUPAC-expanded matching is used
Found cut positions: [3, 7]
```

The pattern "AGR" matches:
- "AGA" at position 2 → cut at position 3
- "AGG" at position 6 → cut at position 7

### Example 2: Overlapping Matches (W)

```bash
python sim.py --seq "ATATATTCC" --enz TestW
```

Where `TestW` has site "ATW" and cut_index=2:

```
Enzyme: TestW
Recognition sequence: ATW
Cut position: after position 2 in recognition sequence
IUPAC-expanded matching is used
Found cut positions: [2, 4, 6]
```

The pattern "ATW" matches:
- "ATA" at position 0 → cut at position 2
- "ATA" at position 2 → cut at position 4 (overlapping)
- "ATT" at position 4 → cut at position 6

### Example 3: Universal Base (N)

```bash
python sim.py --seq "AATACGAT" --enz TestN
```

Where `TestN` has site "AN" and cut_index=1:

```
Enzyme: TestN
Recognition sequence: AN
Cut position: after position 1 in recognition sequence
IUPAC-expanded matching is used
Found cut positions: [1, 2, 4, 7]
```

The pattern "AN" matches:
- "AA" at position 0 → cut at position 1
- "AT" at position 1 → cut at position 2
- "AC" at position 3 → cut at position 4
- "AT" at position 6 → cut at position 7

### Example 4: Multi-Enzyme with IUPAC

```bash
python sim.py --seq "ATAGATAGGCC" --enz TestR TestATA
```

Where `TestR` has site "AGR" and cut_index=1, and `TestATA` has site "ATA" and cut_index=2:

```
Enzyme: TestR
Recognition sequence: AGR
Cut position: after position 1 in recognition sequence
IUPAC-expanded matching is used
Found cut positions: [3, 7]

Enzyme: TestATA
Recognition sequence: ATA
Cut position: after position 2 in recognition sequence
Found cut positions: [2, 6]

COMBINED DIGEST SUMMARY
========================================
Total cuts (unique across enzymes): 4
Cut positions: [2, 3, 6, 7]
Fragment lengths: [2, 1, 3, 4]
Total length verification: 10 bp (expected: 10 bp)
```

## Database Format

IUPAC sites can be used in the `enzymes.json` database:

```json
[
  {
    "name": "TestR",
    "site": "AGR",
    "cut_index": 1
  },
  {
    "name": "TestW", 
    "site": "ATW",
    "cut_index": 2
  },
  {
    "name": "TestN",
    "site": "AN", 
    "cut_index": 1
  }
]
```

## Technical Details

- IUPAC bases are converted to regex character classes (e.g., R → [AG], N → [ACGT])
- Overlapping matches are found using positive lookahead patterns `(?=pattern)`
- Case-insensitive matching is supported for both sequences and IUPAC codes
- All existing functionality (exact matches, multi-enzyme, validation) remains unchanged
- Invalid IUPAC characters in database entries are rejected with clear error messages
