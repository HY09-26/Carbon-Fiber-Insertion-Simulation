"""
config.py
=========
Every setting shared by the notebooks and UI.py, in one place.

The animal list and the probe geometry live here rather than in each
notebook, so every entry point is guaranteed to use the same animals and the
same probes. Edit this file to add a probe or change the insertion protocol.

Paths are resolved relative to this file, so nothing depends on the directory
a notebook is launched from.
"""

import os

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(ROOT, "mouse_data")        # raw tif stacks, not in the repo
SEG_CACHE_DIR = os.path.join(ROOT, "seg_cache")    # cached segmentations, rebuilt on demand
RESULTS_DIR = os.path.join(ROOT, "results")

# Reslice 2 and 3 are excluded: their raw data is too noisy.
ANIMALS = [0, 1, 4, 5, 6, 7, 8, 9]


def tiff_path(animal):
    """Path to one animal's binary vasculature stack."""
    return os.path.join(DATA_DIR, f"Reslice of {animal}.tif")


# Insertion protocol. All lengths are in voxels, which equal micrometres here.
Y_CENTER = 50          # position along the 100-voxel-thick axis: the slab centre
POS_NUM = 100          # insertion sites per animal, evenly spaced along x
DEPTH_LIMIT = 1000     # insertion depth

# Probe geometry, in micrometres (diameters). Each probe is a shank (nearer the
# cortical surface) followed by a tip (deepest), each tapering linearly from
# its base to its top diameter. Keys are used in output filenames.
ELECTRODES = {
    "CF": dict(label="Carbon Fiber 10 um",
               shank_length=840,     shank_base_d=8.4,  shank_top_d=8.4,
               tip_length=160,       tip_base_d=6.8,    tip_top_d=0.0),
    "FMA": dict(label="Microprobes FMA",
                shank_length=971.91, shank_base_d=25.0, shank_top_d=25.0,
                tip_length=28.09,    tip_base_d=25.0,   tip_top_d=0.0),
    "Shuttle": dict(label="Shuttle 25 um",
                    shank_length=990, shank_base_d=25.0, shank_top_d=25.0,
                    tip_length=10,    tip_base_d=25.0,   tip_top_d=25.0),
    "UEA": dict(label="Blackrock UEA",
                shank_length=950,    shank_base_d=90.0, shank_top_d=28.0,
                tip_length=50,       tip_base_d=28.0,   tip_top_d=3.0),
}

# Cylinder control sweep (Bleeding_per_diameter.ipynb): 4, 8, ..., 80 um.
CYLINDER_DIAMETERS = list(range(4, 81, 4))

# 3-D figure. Rendered at a single site; illustrative only, not used in any number.
FIG_X_CENTER = 2200
FIG_CROP_X = 80        # half-width of the rendered slab along x
FIG_CROP_Y = 50        # half-width along y; 50 fills the slab
FIG_FORMAT = "tif"
