#!/usr/bin/env python3
"""
Gel Ladders Module
Defines standard DNA ladder presets for agarose gel simulation.
"""

# Standard ladder presets with fragment sizes in base pairs
LADDER_PRESETS = {
    "100bp": [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000],
    "1kb": [250, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000],
    "broad": [100, 200, 500, 1000, 1500, 2000, 3000, 4000, 5000, 7000, 10000, 12000],
}


def get_ladder(name: str) -> list:
    """
    Get ladder fragment sizes by name.
    
    Args:
        name: Name of the ladder preset ("100bp", "1kb", "broad")
        
    Returns:
        List of fragment sizes in base pairs
        
    Raises:
        ValueError: If ladder name is not recognized
    """
    name_lower = name.lower()
    if name_lower not in LADDER_PRESETS:
        available = ", ".join(LADDER_PRESETS.keys())
        raise ValueError(f"Unknown ladder '{name}'. Available ladders: {available}")
    
    return LADDER_PRESETS[name_lower]


def get_available_ladders() -> list:
    """
    Get list of available ladder names.
    
    Returns:
        List of ladder names
    """
    return list(LADDER_PRESETS.keys())

