#!/usr/bin/env python3
import csv
import json
import re
import sys
import os
from typing import Tuple, List, Dict, Any

IUPAC_ALLOWED = set("ACGTRYSWKMBDHVN")

def normalize_recognition(seq_raw: str) -> str:
    """Keep only IUPAC letters, uppercase (handles 5'-...-3' noise)."""
    if seq_raw is None:
        return ""
    s = seq_raw.strip().upper()
    letters = re.findall(r"[A-Z]", s)
    return "".join(ch for ch in letters if ch in IUPAC_ALLOWED)

def parse_cut_site(cut_raw: str) -> Tuple[str, int]:
    """
    Parse Cut Site like:
      C^GGCCG, AT^CGAT, ^CCWGG, CCWGG^
    If there are alternatives: 'Y^GGCCR OR CGGCCR-C' -> take the first option.
    Return (motif_without_caret, cut_index) where cut_index == #letters before ^.
    """
    if not cut_raw:
        raise ValueError("Cut Site is empty")
    s = cut_raw.strip().upper()
    # Keep only the first alternative if 'OR' is present
    s = re.split(r"\bOR\b", s, maxsplit=1)[0].strip()
    # Remove everything except letters and caret
    cleaned = re.sub(r"[^A-Z^]", "", s)
    if cleaned.count("^") != 1:
        raise ValueError(f"Cut Site must contain exactly one '^': {cut_raw!r}")
    left, right = cleaned.split("^")
    site = left + right
    bad = [ch for ch in site if ch not in IUPAC_ALLOWED]
    if bad:
        raise ValueError(f"Cut Site contains invalid letters {bad} in {cut_raw!r}")
    return site, len(left)

def normalize_overhang_type(raw: str) -> str:
    """
    Normalize overhang type from CSV to standard format.
    
    Args:
        raw: Raw overhang type string from CSV
        
    Returns:
        Normalized overhang type: "5' overhang", "3' overhang", "Blunt", or "Unknown"
    """
    if not raw:
        return "Unknown"
    s = raw.strip().lower().replace("'", "'")
    if "blunt" in s:
        return "Blunt"
    if "5'" in s:
        return "5' overhang"
    if "3'" in s:
        return "3' overhang"
    return "Unknown"

def main():
    if len(sys.argv) < 2:
        print("Usage: python convert_enzyme_table.py <enzymes.csv> [--out enzymes.json]")
        sys.exit(1)

    infile = sys.argv[1]
    out_path = None
    if len(sys.argv) >= 4 and sys.argv[2] == "--out":
        out_path = sys.argv[3]

    if not os.path.exists(infile):
        print(f"ERROR: Could not find CSV file: {infile!r}. If it has spaces, wrap in quotes.")
        sys.exit(2)

    raw_entries: List[Dict[str, Any]] = []
    skipped_rows: List[int] = []

    with open(infile, newline="", encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader, start=2):
            name = (row.get("Enzyme") or "").strip()
            recog = (row.get("Recognition Sequence") or "").strip()
            cut = (row.get("Cut Site") or "").strip()
            overhang_raw = (row.get("Overhang Type") or "").strip()

            if not name and not recog and not cut:
                continue
            if not (name and cut):
                skipped_rows.append(i)
                continue

            recog_seq = normalize_recognition(recog)
            try:
                cut_seq, cut_idx = parse_cut_site(cut)
            except Exception as e:
                print(f"WARNING row {i} {name}: {e}")
                skipped_rows.append(i)
                continue

            # Prefer the cut-site motif if it differs (often includes Ns/spacers)
            site = cut_seq if (not recog_seq or recog_seq != cut_seq) else recog_seq

            # Validate range
            if not (0 <= cut_idx <= len(site)):
                print(f"WARNING row {i} {name}: cut_index {cut_idx} out of range for site {site}")
                skipped_rows.append(i)
                continue

            # Normalize overhang type
            overhang_type = normalize_overhang_type(overhang_raw)
            if overhang_type == "Unknown" and overhang_raw:
                print(f"WARNING row {i} {name}: Unknown overhang type '{overhang_raw}', setting to 'Unknown'")

            raw_entries.append({"name": name, "site": site, "cut_index": cut_idx, "overhang_type": overhang_type})

    # ---- Deduplicate entries (same name + same site + same cut_index + same overhang_type) ----
    seen = set()
    merged: List[Dict[str, Any]] = []
    conflicts: List[str] = []
    
    for entry in raw_entries:
        key = (entry["name"], entry["site"], entry["cut_index"], entry["overhang_type"])
        if key in seen:
            continue
        
        # Check for conflicts (same name + site but different overhang_type)
        conflict_key = (entry["name"], entry["site"])
        existing_overhangs = [e["overhang_type"] for e in merged if (e["name"], e["site"]) == conflict_key]
        if existing_overhangs and entry["overhang_type"] != existing_overhangs[0]:
            conflicts.append(f"{entry['name']} ({entry['site']}): {existing_overhangs[0]} vs {entry['overhang_type']}")
            # Prefer non-Unknown values
            if entry["overhang_type"] != "Unknown" and existing_overhangs[0] == "Unknown":
                # Replace the Unknown entry with the new one
                merged = [e for e in merged if (e["name"], e["site"]) != conflict_key]
            else:
                # Keep existing entry, skip this one
                continue
        
        seen.add(key)
        merged.append({
            "name": entry["name"],
            "site": entry["site"], 
            "cut_index": entry["cut_index"],
            "overhang_type": entry["overhang_type"]
        })

    if out_path is None:
        out_path = os.path.join(os.getcwd(), "enzymes.json")

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(merged, f, indent=2)

    print("\n✅ Conversion complete!")
    print(f"→ Wrote {len(merged)} enzymes to: {out_path}")
    if skipped_rows:
        print(f"⚠️ Skipped {len(skipped_rows)} row(s). Examples: {skipped_rows[:8]}")
    if conflicts:
        print("\n⚠️ Overhang type conflicts found:")
        for conflict in conflicts[:10]:
            print(f"  • {conflict}")
        if len(conflicts) > 10:
            print(f"  ... and {len(conflicts)-10} more conflicts.")

if __name__ == "__main__":
    main()

