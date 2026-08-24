"""
Area_UI.py
==========
DEPRECATED - kept only so that older imports keep working.

This file used to hold its own copy of the volume-metric code.  That copy was
identical to `Volume_bleeding.py` apart from one progress `print`, which meant
the same algorithm lived in two places and could drift apart.  It now simply
re-exports the canonical implementations.

Use `Volume_bleeding` directly in new code:

    from Volume_bleeding import (create_cone_mask, calculate_cone_radius,
                                 simulate_cone_insertion,
                                 find_first_black_pixel_slice,
                                 process_cone_positions)

The only observable difference from the historical version is that
`process_cone_positions` now prints "Counting the volume of intersected
vessels..." while it runs.  Nothing in the pipeline reads that output.
"""

from Volume_bleeding import (
    create_cone_mask,
    calculate_cone_radius,
    simulate_cone_insertion,
    find_first_black_pixel_slice,
    process_cone_positions,
)

__all__ = [
    "create_cone_mask",
    "calculate_cone_radius",
    "simulate_cone_insertion",
    "find_first_black_pixel_slice",
    "process_cone_positions",
]
